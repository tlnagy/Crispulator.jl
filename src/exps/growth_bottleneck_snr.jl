function growth_bottleneck_snr(filepath; debug=false, quiet=false)
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
            :representation => [100],
            :bottleneck_representation => [100],
            :seq_depth => [100],
            :num_bottlenecks => [1, 10]
        )
        num_runs = 1
    end

    runs = grouped_param_space(GrowthScreen(), parameters, [:num_bottlenecks], num_runs)

    max_phenotype_dists = Dict{Symbol, Tuple{Float64, Sampleable}}(
        :inactive => (0.85, Delta(0.0)),
        :negcontrol => (0.05, Delta(0.0)),
        :decreasing => (0.1, TruncatedNormal(-0.55, 0.2, -1, -0.1))
    )

    lib = Library(max_phenotype_dists, CRISPRi())
    before = time()
    results = pmap(args -> run_exp(args[1], lib, compute_snr; run_idx=args[2]), runs)
    (!quiet) && println("$(time() - before) seconds")
    results = DataFrame(permutedims(hcat(results...), [2, 1]))
    results[:crisprtype] = "CRISPRi"
    lib = Library(max_phenotype_dists, CRISPRn())
    before = time()
    results2 = pmap(args -> run_exp(args[1], lib, compute_snr; run_idx=args[2]), runs)
    (!quiet) && println("$(time() - before) seconds")
    results2 = DataFrame(permutedims(hcat(results2...), [2, 1]))
    results2[:crisprtype] = "CRISPRn"
    results = vcat(results, results2)

    hierarchy = reshape([:snr, :signal, :noise], (3, 1))
    new_names = [[:technique, :score]...; fieldnames(GrowthScreen)...; :run_idx; :crisprtype;]

    results = construct_hierarchical_label(hierarchy, results, new_names)
    CSV.write(filepath, results)
end
