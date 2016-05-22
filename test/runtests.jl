using Base.Test

include("../src/load.jl")

println("Running tests:")
filenames = ["kdrelationships.jl"]
for filename in filenames
    include(filename)
    println("\t\033[1m\033[32mPASSED\033[0m: $(filename)")
end

println("All tests pass")
