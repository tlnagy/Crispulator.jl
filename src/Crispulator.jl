module Crispulator

using StatsBase
using Distributions
using DataFrames
using HypothesisTests
using IterTools
using DocStringExtensions
using Combinatorics
using DataStructures

include("common.jl")
include("utils.jl")
include("library.jl")
include("transfection.jl")
include("selection.jl")
include("sequencing.jl")
include("processing.jl")
include("designs.jl")

export FacsScreen, GrowthScreen, CRISPRi, CRISPRn, Delta, Library,
       construct_library, transfect, select, counts_to_freqs, 
       sequencing, differences_between_bins, venn, auprc, auroc, 
       ScreenSetup

end