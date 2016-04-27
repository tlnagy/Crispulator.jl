using DataFrames
using Distributions
using Gadfly
using HypothesisTests

include("utils.jl")
include("common.jl")
include("library.jl")
include("transfection.jl")
include("selection.jl")
include("sequencing.jl")
include("analysis.jl")

const N = 500 # number of target genes
const coverage = 5 # number of guides per gene
const representation = 1000 # Number of cells with each guide

guides, guide_freqs = construct_library(N, coverage)

cell_count = N*coverage*representation
const moi = 0.25
guide_freqs_dist = Categorical(guide_freqs)

cells = transfect(guide_freqs_dist, cell_count, moi)

bin_info = Dict(:bin1 => (0.0, 1/3), :bin2 => (2/3, 1.0))
bin_cells = facs_sort(cells, guides, bin_info, 0.5)

freqs = counts_to_freqs([bin_cells[binname] for binname in keys(bin_cells)]...)
raw_data = sequencing([10^7, 10^7], guides, freqs...)

analyze(raw_data, gen_plots=true)
