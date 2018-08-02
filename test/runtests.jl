fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"
quiet = length(ARGS) > 0 && ARGS[1] == "-q"
errorfound = false

load_file = joinpath("src", "simulation", "load.jl")
include(normpath(joinpath(Base.source_dir(),"..",load_file)))

using Base.Test
using DataStructures
using ColorBrewer
using Gadfly
using Compat
using CSV
import Compat: UTF8String, ASCIIString, view

println("Running tests:")
filenames = ["kdrelationships.jl",
             "qualitymetrics.jl",
             "diffcrisprtransfection.jl",
             "processing.jl",
             "growth.jl",
             "selectionmethods.jl",
             "cmdline.jl"]

for filename in filenames
    try
        include(filename)
        println("\t\033[1m\033[32mPASSED\033[0m:  $(filename)")
    catch e
        errorfound = true
        println("\t\033[1m\033[31mFAILED\033[0m:  $(filename)")
        if fatalerrors
            rethrow(e)
        elseif !quiet
            showerror(STDOUT, e, backtrace())
            println()
        end
    end
end

analyses_path = normpath(joinpath(Base.source_dir(),"..",joinpath("src", "exps")))
analyses = readdir(analyses_path)
analyses = analyses[find(x -> x != "common.jl", analyses)]
# load common file first
include(joinpath(analyses_path, "common.jl"))
skipped = 0

for analysis in analyses
    try
        # skip testing with Gadfly on 32-bit machines, this is an
        # upstream problem
        if analysis == "gen_plots.jl" && Sys.WORD_SIZE != 64
            println("\t\033[1m\033[33mSKIPPED\033[0m: $(analysis)" *
            " # plotting library broken on 32bit machines")
            skipped+=1
            continue
        end
        include(joinpath(analyses_path, analysis))
        tempfile = tempname()
        func_name = Symbol(splitext(analysis)[1])
        getfield(Main, func_name)(tempfile, debug=true, quiet=true)
        println("\t\033[1m\033[32mPASSED\033[0m:  $(analysis)")
    catch e
        errorfound = true
        println("\t\033[1m\033[31mFAILED\033[0m:  $(analysis)")
        if fatalerrors
            rethrow(e)
        elseif !quiet
            showerror(STDOUT, e, backtrace())
            println()
        end
    end
end

if errorfound
    throw("Tests failed")
else
    println("All tests pass, $skipped test$(skipped == 1 ? "" : "s") skipped")
end
