function growth_bottlenecks(filepath; debug=false, quiet=false)
    if !debug
        parameters = Dict{Symbol, Vector}(
            :num_genes => repeat([5000], inner=[3]),
            :representation => [10, 100, 1000],
            :bottleneck_representation => [10, 100, 1000],
            :seq_depth => [100, 100, 1000],
            :num_bottlenecks => collect(1:20),
            :noise => collect(range(0.001, 0.1, length = 5))
        )
        num_runs = 25
    else
        parameters = Dict{Symbol, Vector}(
            :representation => [100],
            :bottleneck_representation => [100],
            :seq_depth => [100],
            :num_bottlenecks => [5, 4],
            :noise => [0.01, 0.5]
        )
        num_runs = 1
    end

    runs = grouped_param_space(GrowthScreen(), parameters, [:num_bottlenecks, :noise], num_runs)

    methods = [venn, auprc]
    measures = [:inc, :dec, :incdec]
    genetypes = [:sigmoidal, :linear, :all]
    test_method_wrapper = (bc_counts, genes) -> test_methods(genes, methods, measures, genetypes)

    max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
        :inactive => (0.849, Delta(0.0)),
        :negcontrol => (0.05, Delta(0.0)),
        :increasing => (0.001, truncated(Normal(0.1, 0.1), 0.025, 1)),
        :decreasing => (0.1, truncated(Normal(-0.55, 0.2), -1, -0.1))
    )

    lib = Library(max_phenotype_dists, CRISPRi())
    p = Progress(length(runs); showspeed = true, enabled = !quiet && !is_logging(stdout))
    results = progress_pmap(args -> Crispulator.run_exp(args[1], lib, test_method_wrapper; run_idx=args[2]), runs; progress = p)
    results = DataFrame(permutedims(hcat(results...), [2, 1]), :auto)
    results[!, :crisprtype] .= "CRISPRi"

    lib = Library(max_phenotype_dists, CRISPRn())
    p = Progress(length(runs); showspeed = true, enabled = !quiet && !is_logging(stdout))
    results2 = progress_pmap(args -> Crispulator.run_exp(args[1], lib, test_method_wrapper; run_idx=args[2]), runs; progress = p)
    results2 = DataFrame(permutedims(hcat(results2...), [2, 1]), :auto)
    results2[!, :crisprtype] .= "CRISPRn"
    results = vcat(results, results2)

    hierarchy = vcat([hcat(item...) for item in Iterators.product(map(Symbol, methods), measures, genetypes)]...)
    new_names = [[:method, :measure, :genetype, :score]...; fieldnames(GrowthScreen)...; :run; :crisprtype]
    results = construct_hierarchical_label(hierarchy, results, new_names)
    CSV.write(filepath, results)
end
