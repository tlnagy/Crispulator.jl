fatalerrors = length(ARGS) > 0 && ARGS[1] == "-f"
quiet = length(ARGS) > 0 && ARGS[1] == "-q"
errorfound = false

using Base.Test
(!isdir(Pkg.dir("DataStructures"))) && Pkg.add("DataStructures")
using DataStructures

include("../src/load.jl")

println("Running tests:")
filenames = ["kdrelationships.jl",
             "qualitymetrics.jl",
             "diffcrisprtransfection.jl",
             "growth.jl",
             "selectionmethods.jl"]
for filename in filenames
    try
        include(filename)
        println("\t\033[1m\033[32mPASSED\033[0m: $(filename)")
    catch e
        errorfound = true
        println("\t\033[1m\033[31mFAILED\033[0m: $(filename)")
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
    println("All tests pass")
end
