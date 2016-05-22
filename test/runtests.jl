using Base.Test
using DataStructures

include("../src/load.jl")

println("Running tests:")
filenames = ["kdrelationships.jl", "diffcrisprtransfection.jl"]
for filename in filenames
    include(filename)
    println("\t\033[1m\033[32mPASSED\033[0m: $(filename)")
end

println("All tests pass")
