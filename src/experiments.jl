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

function grouped_param_space{T <: ScreenSetup}(::T, parameters::Dict{Symbol, Vector}, dist::Symbol, num_runs::Int)
    fields = collect(keys(parameters))
    n_fields = length(fields)
    deleteat!(fields, findin(fields, [dist]))
    runs = []
    grouped_params = zip([parameters[field] for field in fields]...)
    push!(fields, dist)
    for vals in Iterators.product(grouped_params, parameters[dist])
        vals = [vals[1]..., vals[2]]
        for run in 1:num_runs
            setup = FacsScreen()
            for idx in 1:n_fields
                setfield!(setup, fields[idx], vals[idx])
            end
            push!(runs, (setup, run))
        end
    end
    runs
end

function facs_binning(filepath)
    parameters = Dict{Symbol, Vector}(
        :representation => [100, 1000, 1000],
        :bottleneck_representation => [100, 1000, 1000],
        :seq_depth => [100, 1000, 1000],
        :σ => [1.0, 1.0, 0.5],
        :bin_info => [Dict(:bin1 => (0.0, p), :bin2 => (1.0-p, 1.0)) for p in
        [linspace(0.5, 0.1, 30); 0.1.^logspace(0.1, 0.3, 3)]]
    )
    num_runs = 100

    runs = grouped_param_space(FacsScreen(), parameters, :bin_info, num_runs)

    test_methods = (genes, methods, measures, genetypes) -> begin
        local results = []
        for (method, measure, genetype) in Iterators.product(methods, measures, genetypes)
            if genetype != :all
                subgene = genes[genes[:behavior] .== genetype, :]
            else
                subgene = genes
            end
            local result = 0.0
            if measure == :incdec
                result = method(abs(subgene[:pvalmeanprod]), subgene[:class], Set([:increasing, :decreasing]))
            elseif measure == :dec
                result = method(subgene[:pvalmeanprod], subgene[:class], Set([:decreasing]), rev=false)
            else
                result = method(subgene[:pvalmeanprod], subgene[:class], Set([:increasing]))
            end
            (typeof(result) <: Tuple) && (result = result[1])
            push!(results, result)
        end
        (results...)
    end

    methods = [venn, auprc]
    measures = [:inc, :dec, :incdec]
    genetypes = [:sigmoidal, :linear, :all]
    test_method_wrapper = genes -> test_methods(genes, methods, measures, genetypes)

    results = @time pmap(args -> run_exp(args[1], Library(CRISPRi()), test_method_wrapper; run_idx=args[2]), runs)
    results = DataFrame(hcat(results...)')
    results[:crisprtype] = "CRISPRi"
    results2 = @time pmap(args -> run_exp(args[1], Library(CRISPRKO()), test_method_wrapper; run_idx=args[2]), runs)
    results2 = DataFrame(hcat(results2...)')
    results2[:crisprtype] = "CRISPRKO"
    results = vcat(results, results2)

    col_names = map(x->symbol(x[1],"_",x[2], "_", x[3]), collect(Iterators.product(map(symbol, methods), measures, genetypes)))
    names!(results, [col_names...; fieldnames(FacsScreen)...; :run; :crisprtype])
    results[:bin_info] = Float64[el[:bin1][2] for el in results[:bin_info]]
    writetable(filepath, results)
end

function scan_perf_of_diff_cat_in_growth(filepath)

    parameters = Dict{Symbol, Vector}(
        :representation => map(x->round(Int64, x), logspace(0, 3, 2)),
        :bottleneck_representation => map(x->round(Int64, x),  logspace(0,3,2)),
        :num_bottlenecks => collect(2:10:20)
    )
    runs = build_parameter_space(GrowthScreen(), parameters, 1)

    max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
        :inactive => (0.83, Delta(0.0)),
        :negcontrol => (0.05, Delta(0.0)),
        :increasing => (0.02, TruncatedNormal(0.1, 0.1, 0.025, 1)),
        :decreasing => (0.1, TruncatedNormal(-0.55, 0.2, -1, -0.1))
    )
    lib = Library(max_phenotype_dists, CRISPRi())

    # computes the aurocs of increasing and decreasing genes separately
    get_aurocs = genes -> begin
        d = auroc(genes[:pvalmeanprod], genes[:class], Set([:decreasing]))
        i = auroc(genes[:pvalmeanprod], genes[:class], Set([:increasing]))
        (d[1], i[1])
    end

    results = @time pmap(args -> run_exp(args[1], lib, get_aurocs; run_idx=args[2]), runs)
    results = DataFrame(hcat(results...)')

    names!(results, [:decreasing; :increasing; fieldnames(GrowthScreen)...; :run])
    writetable(filepath, results)
end

function scan_representation_space(filepath)
    parameters = Dict{Symbol, Vector}(
        :representation => map(x->round(Int64, x), logspace(0, 3, 20)),
        :bottleneck_representation => map(x->round(Int64, x),  logspace(0,3,20)),
        :seq_depth => [100, 1000, 10000]
    )
    num_runs = 25

    const overlap = intersect(fieldnames(FacsScreen), fieldnames(GrowthScreen))
    # custom function for handling both growth and a facs screen in the
    # same relational datastructure
    flatten_overlap = (screen) -> begin
        local results = Any[typeof(screen)]
        for name in overlap
            push!(results, getfield(screen, name))
        end
        results
    end

    # computes the auprcs of increasing and decreasing genes separately
    get_auprcs = genes -> begin
        a = auprc(abs(genes[:pvalmeanprod]), genes[:class], Set([:increasing, :decreasing]))
        d = auprc(genes[:pvalmeanprod], genes[:class], Set([:decreasing]), rev=false)
        i = auprc(genes[:pvalmeanprod], genes[:class], Set([:increasing]))
        (a[1], d[1], i[1])
    end

    results = []
    for screentype in [FacsScreen(), GrowthScreen()]
        for crisprtype in [CRISPRi(), CRISPRKO()]
            if typeof(screentype) == GrowthScreen
                max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
                    :inactive => (0.83, Delta(0.0)),
                    :negcontrol => (0.05, Delta(0.0)),
                    :increasing => (0.02, TruncatedNormal(0.1, 0.1, 0.025, 1)),
                    :decreasing => (0.1, TruncatedNormal(-0.55, 0.2, -1, -0.1))
                )
                lib = Library(max_phenotype_dists, crisprtype)
            else
                lib = Library(crisprtype)
            end

            runs = build_parameter_space(screentype, parameters, num_runs)

            result = @time pmap(args -> run_exp(args[1], lib, get_auprcs; run_idx=args[2], flatten_func=flatten_overlap), runs)
            result = DataFrame(hcat(result...)')
            result[:crisprtype] = typeof(crisprtype)
            push!(results, result)
        end
    end
    results = vcat(results...)
    names!(results, [:all; :inc; :dec; :screen; overlap...; :run; :crisprtype])
    writetable(filepath, results)
end
