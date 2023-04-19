using Documenter
using Crispulator

makedocs(
    modules = [Crispulator],
    clean = false,
    format = Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        assets =[
                asset("https://analytics.tamasnagy.com/js/script.js", class=:js, attributes=Dict(Symbol("data-domain") => "tamasnagy.com", :defer => ""))
            ],
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
