using Crispulator

using Test
using DataFrames
using DataStructures
using Gadfly
using CSV
using Distributions
using StatsBase
using Random

filenames = ["kdrelationships.jl",
             "qualitymetrics.jl",
             "diffcrisprtransfection.jl",
             "processing.jl",
             "growth.jl",
             "selectionmethods.jl",
             "cmdline.jl"]

for filename in filenames
    @testset "$filename" begin
        include(filename)
    end
end

analyses_path = joinpath(@__DIR__, "..", "exps")
analyses = readdir(analyses_path)
analyses = analyses[findall(x -> x != "common.jl", analyses)]
# load common file first
include(joinpath(analyses_path, "common.jl"))
skipped = 0

for analysis in analyses
    # skip testing with Gadfly on 32-bit machines, this is an
    # upstream problem
    if analysis == "gen_plots.jl" && Sys.WORD_SIZE != 64
        println("\t\033[1m\033[33mSKIPPED\033[0m: $(analysis)" *
        " # plotting library broken on 32bit machines")
        continue
    end
    include(joinpath(analyses_path, analysis))
    tempfile = tempname()
    func_name = Symbol(splitext(analysis)[1])
    @testset "$analysis" begin
        getfield(Main, func_name)(tempfile, debug=true, quiet=true)
    end
end