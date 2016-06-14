print("Loading simulation framework...")

packages = [:DataFrames,
            :Distributions,
            :HypothesisTests,
            :StatsBase,
            :Iterators]
for package in packages
    !(isdir(Pkg.dir(string(package)))) && Pkg.add(string(package))
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

println("Done.")
