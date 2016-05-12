packages = [:DataFrames,
            :Distributions,
            :HypothesisTests,
            :StatsBase,
            :Iterators]
all_workers = packages[1:4]
for package in packages
    try
        Pkg.installed(string(package))
        eval(:(using $package))
        if package in all_workers
            eval(:(@everywhere using $package))
        end
    catch
        Pkg.add(string(package))
    end
end

# load all simulations files on all workers
@everywhere function include_all()
    filenames = ["common.jl", "utils.jl", "library.jl", "transfection.jl",
                 "selection.jl", "sequencing.jl", "processing.jl",
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
