using DataFrames
using Distributions
using HypothesisTests
using StatsBase
@everywhere using DataFrames
@everywhere using Distributions
@everywhere using HypothesisTests
@everywhere using StatsBase

@everywhere function include_all()
    filenames = ["common.jl", "utils.jl", "library.jl", "transfection.jl",
                 "selection.jl", "sequencing.jl", "analysis.jl"]
    for filename in filenames
        include(filename)
    end
end
@everywhere include_all()

@everywhere function run_exp()
    num_genes = 500 # number of target genes
    coverage = 5 # number of guides per gene
    representation = 100 # Number of cells with each guide
    moi = 0.25 # multiplicity of infection
    σ = 1.0 # std dev expected for cells during facs sorting (in phenotype units)
    bin_info = Dict(:bin1 => (0.0, 1/3), :bin2 => (2/3, 1.0))

    guides, guide_freqs = construct_library(num_genes, coverage)

    guide_count = num_genes * coverage
    cell_count = guide_count*representation

    guide_freqs_dist = Categorical(guide_freqs)

    cells = transfect(guides, guide_freqs_dist, cell_count, moi)

    bin_cells = facs_sort(cells, guides, bin_info, σ)

    freqs = counts_to_freqs(bin_cells, guide_count)
    raw_data = sequencing(Dict(:bin1=>10^7,:bin2=>10^7), guides, freqs)

    auroc = analyze(raw_data, gen_plots=false)

    auroc
end

@everywhere println(run_exp())
