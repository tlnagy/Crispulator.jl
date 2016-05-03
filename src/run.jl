using DataFrames
using Distributions
using Gadfly
using HypothesisTests

include("common.jl")
include("utils.jl")
include("library.jl")
include("transfection.jl")
include("selection.jl")
include("sequencing.jl")
include("analysis.jl")

N = 500 # number of target genes
coverage = 5 # number of guides per gene
representation = 100 # Number of cells with each guide
moi = 0.25 # multiplicity of infection
σ = 1.0 # std dev expected for cells during facs sorting (in phenotype units)
bin_info = Dict(:bin1 => (0.0, 1/3), :bin2 => (2/3, 1.0))

function run_facs_crispri(num_genes::Int64,
                          coverage::Int64,
                          representation::Int64,
                          moi::Float64,
                          σ::Float64,
                          bin_info::Dict{Symbol, Tuple{Float64, Float64}})

    guides, guide_freqs = construct_library(num_genes, coverage)

    cell_count = num_genes*coverage*representation
    guide_freqs_dist = Categorical(guide_freqs)

    cells = transfect(guides, guide_freqs_dist, cell_count, moi)

    bin_cells = facs_sort(cells, guides, bin_info, σ)

    freqs = counts_to_freqs(bin_cells)
    raw_data = sequencing(Dict(:bin1=>10^7,:bin2=>10^7), guides, freqs)

    auroc = analyze(raw_data, gen_plots=true)

    auroc
end

auroc = run_facs_crispri(N, coverage, representation, moi, σ, bin_info)
