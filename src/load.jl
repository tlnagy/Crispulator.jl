info("Loading simulation framework")

packages = [:StatsBase,
            :Distributions,
            :DataFrames,
            :HypothesisTests,
            :IterTools,
            :DocStringExtensions,
            :Combinatorics,
            :DataStructures]

for package in packages
    eval(:(using $package))
    eval(:(@everywhere using $package))
end

# Load all simulations files on all workers
@everywhere function include_all(curr_dir)
    filenames = ["common.jl", "utils.jl", "library.jl", "transfection.jl",
                 "selection.jl", "sequencing.jl", "processing.jl",
                 "designs.jl"]
    for filename in filenames
        include(joinpath(curr_dir, filename))
    end
end
curr_dir = Base.source_dir()
@eval @everywhere include_all($curr_dir)
