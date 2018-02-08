"""
$(SIGNATURES)

Runs a screen given the parameters specified in `setup` using the
library `lib` and applies the `processing_func` function to the result.
"""
function run_exp(setup::ScreenSetup, lib::Library, processing_func::Function; run_idx=-1, flatten_func::Function=flatten_setup)

    guides, guide_freqs_dist = construct_library(setup, lib)

    cells, cell_phenotypes = transfect(setup, lib, guides, guide_freqs_dist)

    bin_cells = select(setup, cells, cell_phenotypes, guides)

    freqs = counts_to_freqs(bin_cells, length(guides))
    # uniform for now for all bins
    # work around for dict comprehension syntax change in v0.4->v0.5
    # https://github.com/JuliaLang/Compat.jl/issues/231
    seq_depths = Dict{Symbol, Int}()
    for binname in keys(bin_cells)
        seq_depths[binname] = setup.seq_depth
    end
    raw_data = sequencing(seq_depths, guides, freqs)
    bc_counts, genes = differences_between_bins(raw_data)

    [processing_func(bc_counts, genes)...; flatten_func(setup, lib)...; run_idx]
end

flatten_setup(setup::ScreenSetup, ::Library) = as_array(setup)
flatten_both(setup::ScreenSetup, lib::Library) = [as_array(setup)...; as_array(lib)...]
