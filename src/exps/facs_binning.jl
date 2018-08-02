function facs_binning(filepath; debug=false, quiet=false)
    if !debug
        parameters = Dict{Symbol, Vector}(
            :representation => [10, 100, 1000],
            :bottleneck_representation => [10, 100, 1000],
            :seq_depth => [10, 100, 1000],
            :σ => [1.0, 1.0, 0.5],
            :bin_info => [OrderedDict(:bin1 => (0.0, p), :bin2 => (1.0-p, 1.0)) for p in linspace(0.5, 0.025, 30)]
        )
        num_runs = 100
    else
        parameters = Dict{Symbol, Vector}(
            :representation => [100],
            :bottleneck_representation => [100],
            :seq_depth => [100],
            :σ => [1.0],
            :bin_info => [OrderedDict(:bin1 => (0.0, 0.25), :bin2 => (0.75, 1.0))]
        )
        num_runs = 1
    end

    runs = grouped_param_space(FacsScreen(), parameters, [:bin_info], num_runs)

    methods = [venn, auprc]
    measures = [:inc, :dec, :incdec]
    genetypes = [:sigmoidal, :linear, :all]
    test_method_wrapper = (bc_counts, genes) -> test_methods(genes, methods, measures, genetypes)

    before = time()
    results = pmap(args -> run_exp(args[1], Library(CRISPRi()), test_method_wrapper; run_idx=args[2]), runs)
    (!quiet) && println("$(time() - before) seconds")
    results = DataFrame(permutedims(hcat(results...), [2, 1]))
    results[:crisprtype] = "CRISPRi"
    before = time()
    results2 = pmap(args -> run_exp(args[1], Library(CRISPRn()), test_method_wrapper; run_idx=args[2]), runs)
    (!quiet) && println("$(time() - before) seconds")
    results2 = DataFrame(permutedims(hcat(results2...), [2, 1]))
    results2[:crisprtype] = "CRISPRn"
    results = vcat(results, results2)

    hierarchy = vcat([hcat(item...) for item in IterTools.product(map(Symbol, methods), measures, genetypes)]...)
    new_names = [[:method, :measure, :genetype, :score]...; fieldnames(FacsScreen)...; :run; :crisprtype]
    results = construct_hierarchical_label(hierarchy, results, new_names)
    results[:bin_info] = Float64[el[:bin1][2] for el in results[:bin_info]]
    CSV.write(filepath, results)
end
