function facs_binning_snr(filepath; debug=false, quiet=false)
    if !debug
        parameters = Dict{Symbol, Vector}(
            :representation => [10, 100, 1000],
            :bottleneck_representation => [10, 100, 1000],
            :seq_depth => [10, 100, 1000],
            :σ => [1.0, 1.0, 0.5],
            :bin_info => [OrderedDict(:bin1 => (0.0, p), :bin2 => (1.0-p, 1.0)) for p in range(0.5, 0.025, length=30)]
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

    p = Progress(length(runs); showspeed = true, enabled = !quiet && !is_logging(stdout))
    results = progress_map(args -> Crispulator.run_exp(args[1], Library(CRISPRi()), compute_snr; run_idx=args[2]), runs, progress = p)
    results = DataFrame(permutedims(hcat(results...), [2, 1]), :auto)
    results[!, :crisprtype] .= "CRISPRi"

    p = Progress(length(runs); showspeed = true, enabled = !quiet && !is_logging(stdout))
    results2 = progress_pmap(args -> Crispulator.run_exp(args[1], Library(CRISPRn()), compute_snr; run_idx=args[2]), runs, progress = p)
    results2 = DataFrame(permutedims(hcat(results2...), [2, 1]), :auto)
    results2[!, :crisprtype] .= "CRISPRn"
    results = vcat(results, results2)

    hierarchy = reshape([:snr, :signal, :noise], (3, 1))
    new_names = [[:technique, :score]...; fieldnames(FacsScreen)...; :run_idx; :crisprtype;]

    results = construct_hierarchical_label(hierarchy, results, new_names)
    results[!, :bin_info] = Float64[el[:bin1][2] for el in results[!, :bin_info]]
    CSV.write(filepath, results)
end
