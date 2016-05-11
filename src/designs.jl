function setup_screen(setup::ScreenSetup, lib::Library)
    guides, guide_freqs = construct_library(lib, setup.num_genes, setup.coverage)
    guide_freqs_dist = Categorical(guide_freqs)
    guides, guide_freqs_dist
end

function teardown_screen(setup::ScreenSetup,
                         guides::Vector{Barcode},
                         bin_cells::Dict{Symbol,Vector{Int64}},
                         testalg::Bool)

    freqs = counts_to_freqs(bin_cells, length(guides))
    raw_data = sequencing(Dict(:bin1=>setup.seq_depth,:bin2=>setup.seq_depth), guides, freqs)

    analyze(raw_data, gen_plots=false, testmethod=testalg)

end

function run_exp(setup::FacsScreen, lib::Library; run_idx=-1, testalg=false)

    guides, guide_freqs_dist = setup_screen(setup, lib)

    min_perc = minimum([range[2] - range[1] for (binname, range) in setup.bin_info])
    expand_to = round(Int64, setup.num_cells_per_bin/min_perc)

    cells = transfect(guides, guide_freqs_dist, length(guides)*setup.representation, setup.moi, expand_to)

    bin_cells = facs_sort(cells, guides, setup.bin_info, setup.σ)

    [teardown_screen(setup, guides, bin_cells, testalg)...; as_array(setup)...; run_idx]
end

function run_exp(setup::GrowthScreen, lib::Library; run_idx=-1, testalg=false)

    guides, guide_freqs_dist = setup_screen(setup, lib)

    cells, num_doublings = transfect(setup, guides, guide_freqs_dist)

    bin_cells = growth_assay(cells, guides, setup.num_bottlenecks, setup.bottleneck_representation)

    [teardown_screen(setup, guides, bin_cells, testalg)...; as_array(setup)...; run_idx]
end
