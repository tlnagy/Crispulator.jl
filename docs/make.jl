using Documenter

makedocs(
    clean = false,
    format = Documenter.Formats.HTML,
    sitename = "Crispulator.jl",
    pages = Any[
        "Home" => "index.md",
        "FACS vs Growth" => "facs_growth.md"
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
