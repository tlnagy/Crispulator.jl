# Growth screens seem to be highly susceptible to the distribution of
# genes that promote growth versus that hurt it. This experiment tests
# the behavior of growth screens under several different library setups

function growth_sensitivity_library(filepath; debug=false, quiet=false)

    if !debug

        inc_space = range(0.02, 0.1, length=5)
        cent_space = range(0.1, 0.55, length = 5)

        libs = Library[]
        for (inc, cent, crisprtype) in Iterators.product(inc_space, cent_space, [CRISPRi(), CRISPRn()])

            max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
                :inactive => (1-inc-0.1-0.05, Delta(0.0)),
                :negcontrol => (0.05, Delta(0.0)),
                :increasing => (inc, truncated(Normal(cent, 0.1), 0.025, 1)),
                :decreasing => (0.1, truncated(Normal(-0.55, 0.2), -1, -0.1))
            )
            push!(libs, Library(max_phenotype_dists, crisprtype))
        end

        parameters = Dict{Symbol, Vector}(
                :representation => [100],
                :bottleneck_representation => [100],
                :seq_depth => [100],
                :num_bottlenecks => collect(1:20)
            )
        num_runs = 25
    else
        libs = Library[Library(CRISPRi()), Library(CRISPRn())]

        parameters = Dict{Symbol, Vector}(
                :representation => [100],
                :bottleneck_representation => [100],
                :seq_depth => [100],
                :num_bottlenecks => [7, 20]
            )
        num_runs = 1
    end

    runs = grouped_param_space(GrowthScreen(), parameters, libs, [:num_bottlenecks], num_runs);

    function compute_auprc(barcodes, genes)
        i = auprc(genes[!, :pvalmeanprod_bin2_div_bin1], genes[!, :class], Set([:increasing]))[1]
        d = auprc(genes[!, :pvalmeanprod_bin2_div_bin1], genes[!, :class], Set([:decreasing]), rev=false)[1]
        (i, d)
    end

    p = Progress(length(runs); showspeed = true, enabled = !quiet && !is_logging(stdout))
    results = progress_pmap(args -> Crispulator.run_exp(args[1], args[2], compute_auprc;
                   run_idx=args[3], flatten_func=Crispulator.flatten_both), runs, progress = p)
    @info "Constructing CSV"
    data = DataFrame(permutedims(hcat(results...), [2, 1]), :auto)
    hierarchy = reshape([:inc; :dec], (2, 1))
    new_names = [[:measure, :score]...; fieldnames(GrowthScreen)...; Crispulator.array_names(Library)...; :run_idx]

    data = construct_hierarchical_label(hierarchy, data, new_names)

    CSV.write(filepath, data)
end
