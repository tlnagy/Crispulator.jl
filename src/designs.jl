"""
Runs a screen given the parameters specified in `setup` using the
library `lib` and applies the `processing_func` function to the result.
"""
function run_exp(setup::ScreenSetup, lib::Library, processing_func::Function; run_idx=-1)

    guides, guide_freqs_dist = construct_library(setup, lib)

    cells, cell_phenotypes = transfect(setup, lib, guides, guide_freqs_dist)

    bin_cells = select(setup, cells, cell_phenotypes, guides)

    freqs = counts_to_freqs(bin_cells, length(guides))
    # uniform for now for all bins
    seq_depths = Dict{Symbol, Int64}([binname=>setup.seq_depth for binname in keys(bin_cells)])
    raw_data = sequencing(seq_depths, guides, freqs)
    genes = differences_between_bins(raw_data)

    [processing_func(genes)...; as_array(setup)...; run_idx]
end
