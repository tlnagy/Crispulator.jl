function run_exp(setup::ScreenSetup, lib::Library; run_idx=-1, testalg=false)

    guides, guide_freqs = construct_library(lib, setup.num_genes, setup.coverage)

    guide_count = setup.num_genes * setup.coverage
    cell_count = guide_count*setup.representation

    guide_freqs_dist = Categorical(guide_freqs)

    min_perc = minimum([range[2] - range[1] for (binname, range) in setup.bin_info])
    expand_to = round(Int64, setup.num_cells_per_bin/min_perc)

    cells = transfect(guides, guide_freqs_dist, cell_count, setup.moi, expand_to)

    bin_cells = facs_sort(cells, guides, setup.bin_info, setup.σ)

    freqs = counts_to_freqs(bin_cells, guide_count)
    raw_data = sequencing(Dict(:bin1=>setup.seq_depth,:bin2=>setup.seq_depth), guides, freqs)

    auroc = analyze(raw_data, gen_plots=false, testmethod=testalg)

    [auroc...; as_array(setup)...; run_idx]
end
