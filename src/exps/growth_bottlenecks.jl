function main(filepath; debug=false, quiet=false)
    if !debug
        parameters = Dict{Symbol, Vector}(
            :representation => [100, 1000, 1000],
            :bottleneck_representation => [100, 1000, 1000],
            :seq_depth => [100, 1000, 1000],
            :num_bottlenecks => collect(1:20)
        )
        num_runs = 25
    else
        parameters = Dict{Symbol, Vector}(
            :representation => [100],
            :bottleneck_representation => [100],
            :seq_depth => [100],
            :num_bottlenecks => [5]
        )
        num_runs = 1
    end

    runs = grouped_param_space(GrowthScreen(), parameters, :num_bottlenecks, num_runs)

    methods = [venn, auprc]
    measures = [:inc, :dec, :incdec]
    genetypes = [:sigmoidal, :linear, :all]
    test_method_wrapper = (bc_counts, genes) -> test_methods(genes, methods, measures, genetypes)

    before = time()
    results = pmap(args -> run_exp(args[1], Library(CRISPRi()), test_method_wrapper; run_idx=args[2]), runs)
    (!quiet) && println("$(time() - before) seconds")
    results = DataFrame(hcat(results...)')
    results[:crisprtype] = "CRISPRi"
    before = time()
    results2 = pmap(args -> run_exp(args[1], Library(CRISPRKO()), test_method_wrapper; run_idx=args[2]), runs)
    (!quiet) && println("$(time() - before) seconds")
    results2 = DataFrame(hcat(results2...)')
    results2[:crisprtype] = "CRISPRKO"
    results = vcat(results, results2)

    hierarchy = vcat([hcat(item...) for item in Iterators.product(map(symbol, methods), measures, genetypes)]...)
    new_names = [[:method, :measure, :genetype, :score]...; fieldnames(GrowthScreen)...; :run; :crisprtype]
    results = construct_hierarchical_label(hierarchy, results, new_names)
    writetable(filepath, results)
end
