using Documenter

module Simulation
    packages = [:StatsBase,
                :Distributions,
                :DataFrames,
                :HypothesisTests,
                :IterTools,
                :DocStringExtensions]

    for package in packages
        eval(:(using $package))
    end

    filenames = ["common.jl", "utils.jl", "library.jl", "transfection.jl",
                 "selection.jl", "sequencing.jl", "processing.jl",
                 "designs.jl"]
    for filename in filenames
        include(joinpath(Base.source_dir(), "..", "src", "simulation", filename))
    end
end
using Simulation

makedocs(
    modules = [Simulation],
    clean = false,
    format = :html,
    sitename = "Crispulator.jl",
    pages = Any[
        "Home" => "index.md",
        "Simulation Internals" => "internals.md"
    ]
)

deploydocs(
    repo   = "github.com/tlnagy/Crispulator.jl.git",
    julia  = "release",
    osname = "linux",
    deps = nothing,
    make = nothing,
    target = "build",
)
