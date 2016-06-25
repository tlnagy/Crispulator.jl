function main(filepath; debug=false, quiet=false)
    if !debug
        parameters = Dict{Symbol, Vector}(
            :representation => map(x->round(Int64, x), logspace(0, 4, 12)),
            :bottleneck_representation => map(x->round(Int64, x),  logspace(0,4,12)),
            :seq_depth => map(x->round(Int64, x),  logspace(0,4,12)),
            :num_bottlenecks => collect(8:2:20)
        )
        num_runs = 10
    else
        parameters = Dict{Symbol, Vector}(
            :representation => map(x->round(Int64, x), logspace(0, 2, 2)),
            :bottleneck_representation => map(x->round(Int64, x),  logspace(0,2,2)),
            :seq_depth => [10^2],
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
    libs = [Library(max_phenotype_dists, CRISPRi()), Library(max_phenotype_dists, CRISPRKO())]
    runs = build_parameter_space(GrowthScreen(), parameters, libs, num_runs);

    before = time()
    results = pmap(args -> run_exp(args[1], args[2], test_method_wrapper;
                   run_idx=args[3], flatten_func=flatten_both), runs)
    (!quiet) && println("$(time() - before) seconds")

    data = DataFrame(hcat(results...)')
    hierarchy = vcat([hcat(item...) for item in Iterators.product(map(symbol, methods), measures, genetypes)]...)
    new_names = [[:method, :measure, :genetype, :score]...; fieldnames(GrowthScreen)...; array_names(Library)...; :run_idx]

    data = construct_hierarchical_label(hierarchy, data, new_names)
    writetable(filepath, data)
end
