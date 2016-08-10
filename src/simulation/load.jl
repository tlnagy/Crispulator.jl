"""
A pooled genetic screen simulator
"""
module CRISPulator

export ScreenSetup, FacsScreen, GrowthScreen
export Library, CRISPRi, CRISPRn
export KDPhenotypeRelationship
export run_exp
export auroc, venn, auprc, signal, noise

using StatsBase
using Distributions
using DataFrames
using HypothesisTests
using Iterators
using Compat
import Compat: UTF8String, ASCIIString, view

include("common.jl")
include("utils.jl")
include("library.jl")
include("transfection.jl")
include("selection.jl")
include("sequencing.jl")
include("processing.jl")
include("designs.jl")

end
