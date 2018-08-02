function growth_representation(filepath; debug=false, quiet=false)
    if !debug
        parameters = Dict{Symbol, Vector}(
            :representation => [10, 100, 1000],
            :bottleneck_representation => [10, 100, 1000],
            :seq_depth => [10, 100, 1000],
            :num_bottlenecks => collect(1:20)
        )
        num_runs = 25
    else
        parameters = Dict{Symbol, Vector}(
            :representation => [1, 10],
            :bottleneck_representation => [1, 10],
            :seq_depth => [1, 10],
            :num_bottlenecks => [10]
        )
        num_runs = 1
    end

    methods = [venn, auprc]
    measures = [:dec]
    genetypes = [:sigmoidal, :linear, :all]
    test_method_wrapper = (bc_counts, genes) -> test_methods(genes, methods, measures, genetypes)

    max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
        :inactive => (1-0.1-0.05, Delta(0.0)),
        :negcontrol => (0.05, Delta(0.0)),
        :decreasing => (0.1, TruncatedNormal(-0.55, 0.2, -1, -0.1))
    )
    libs = [Library(max_phenotype_dists, CRISPRi()), Library(max_phenotype_dists, CRISPRn())]
    runs = grouped_param_space(GrowthScreen(), parameters, libs, [:num_bottlenecks], num_runs);

    before = time()
    results = pmap(args -> run_exp(args[1], args[2], test_method_wrapper;
                   run_idx=args[3], flatten_func=flatten_both), runs)
    (!quiet) && println("$(time() - before) seconds")

    data = DataFrame(permutedims(hcat(results...), [2, 1]))
    hierarchy = vcat([hcat(item...) for item in IterTools.product(map(Symbol, methods), measures, genetypes)]...)
    new_names = [[:method, :measure, :genetype, :score]...; fieldnames(GrowthScreen)...; array_names(Library)...; :run_idx]

    data = construct_hierarchical_label(hierarchy, data, new_names)
    CSV.write(filepath, data)
end
