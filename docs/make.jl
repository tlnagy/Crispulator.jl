using Documenter

"""
Wrapper for src/simulation code. This is duplicitous since this should really
be inside of src/simulation/load.jl but I had issues with threading and module
loading in the past so this approach should work for now even if it is clunky.
"""
module Simulation
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
    end

    filenames = ["common.jl", "utils.jl", "library.jl", "transfection.jl",
                 "selection.jl", "sequencing.jl", "processing.jl",
                 "designs.jl"]
    for filename in filenames
        include(joinpath(Base.source_dir(), "..", "src", "simulation", filename))
    end

    export FacsScreen, GrowthScreen, Library, Sampleable, Delta, CRISPRi, CRISPRn,
    construct_library, transfect, select, counts_to_freqs, sequencing,
    differences_between_bins, auprc
end
using Simulation

makedocs(
    modules = [Simulation],
    clean = false,
    format = :html,
    sitename = "Crispulator.jl",
    pages = Any[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Custom Simulations" => "custom.md",
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
