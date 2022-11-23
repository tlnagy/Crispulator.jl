using Documenter
using Crispulator

makedocs(
    modules = [Crispulator],
    clean = false,
    format = Documenter.HTML(
        prettyurls = true
    ),
    sitename = "Crispulator.jl",
    pages = [
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "Custom Simulations" => "custom.md",
        "Multiprocessing" => "multiprocessing.md",
        "Crispulator Internals" => "internals.md"
    ]
)

deploydocs(
    repo   = "github.com/tlnagy/Crispulator.jl.git",
    push_preview = true
)
