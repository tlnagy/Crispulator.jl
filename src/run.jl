using DataFrames
using Distributions
using HypothesisTests
using StatsBase
using Iterators

# these packages are used in the simulations and must be loaded on
# all workers
@everywhere using DataFrames
@everywhere using Distributions
@everywhere using HypothesisTests
@everywhere using StatsBase

# load all simulations files on all workers
@everywhere function include_all()
    filenames = ["common.jl", "utils.jl", "library.jl", "transfection.jl",
                 "selection.jl", "sequencing.jl", "analysis.jl",
                 "designs.jl"]
    for filename in filenames
        include(filename)
    end
end
@everywhere include_all()

include("experiments.jl")

function main()
    target_function = symbol(ARGS[1])
    filename = ARGS[2]
    println("Calling $target_function and saving in $filename")
    eval(Expr(:call, target_function, filename))
end

# fire up simulation if run using command line
if !isinteractive()
    # run_wrapper("../data/test_alg.csv")
    main()
end
