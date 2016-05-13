function build_parameter_space{T <: ScreenSetup}(::T, parameters::Dict{Symbol, Vector}, num_runs::Int)
    fields = collect(keys(parameters))
    n_fields = length(fields)
    runs = []
    for vals in Iterators.product([parameters[field] for field in fields]...)
        for run in 1:num_runs
            setup = T()
            for idx in 1:n_fields
                setfield!(setup, fields[idx], vals[idx])
            end
            push!(runs, (setup, run))
        end
    end
    runs
end

function scan_best_methods(filepath)
    parameters = Dict{Symbol, Vector}(
        :representation => map(x->round(Int64, x), logspace(0, 3, 2)),
        :num_cells_per_bin => map(x->round(Int64, x),  2.5*logspace(3,6,2)),
        :seq_depth => map(x->round(Int64, x),  logspace(5,7,2)),
        :σ => [0.5, 1, 2]
    )
    runs = build_parameter_space(FacsScreen(), parameters, 10)

    results = @time pmap(args -> run_exp(args[1], Library(); run_idx=args[2], testalg=true), runs)
    results = DataFrame(hcat(results...)')

    names!(results, [:pval_auroc; :mean_auroc; :pvalmeanprod_auroc ; fieldnames(FacsScreen)...; :run])
    # remove bin_info for now because there isn't a good way to encode
    # that information
    delete!(results, :bin_info)
    writetable(filepath, results)
end

function scan_perf_of_diff_cat_in_growth(filepath)

    parameters = Dict{Symbol, Vector}(
        :representation => map(x->round(Int64, x), logspace(0, 3, 10)),
        :bottleneck_representation => map(x->round(Int64, x),  logspace(0,3,10)),
        :num_bottlenecks => map(x->round(Int64, x), linspace(1, 25, 10)),
        :seq_depth => map(x->round(Int64, x),  logspace(5,7,10))
    )
    runs = build_parameter_space(GrowthScreen(), parameters, 25)

    # computes the aurocs of increasing and decreasing genes separately
    get_aurocs = genes -> begin
        notincreasing = genes[genes[:class] .!= :increasing, :]
        notdecreasing = genes[genes[:class] .!= :decreasing, :]
        d = compute_auroc(Vector(notincreasing[:pvalmeanprod]), Vector(notincreasing[:class]), 100)
        i = compute_auroc(Vector(notdecreasing[:pvalmeanprod]), Vector(notdecreasing[:class]), 100)
        (d, i)
    end

    results = @time pmap(args -> run_exp(args[1], Library(), get_aurocs; run_idx=args[2]), runs)
    results = DataFrame(hcat(results...)')

    names!(results, [:decreasing; :increasing; fieldnames(GrowthScreen)...; :run])
    writetable(filepath, results)
end
