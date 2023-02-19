cd(@__DIR__)

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/apply_style.jl",
    joinpath(@__DIR__, "apply_style.jl")
)
include("apply_style.jl")

using TransitionIndicators

DYNAMICALSYSTEMSBASE_PAGES = [
    "index.md",
]

makedocs(
    modules = [TransitionIndicators],
    format = Documenter.HTML(
        prettyurls = CI,
        assets = [
            asset("https://fonts.googleapis.com/css?family=Montserrat|Source+Code+Pro&display=swap", class=:css),
        ],
    ),
    sitename = "TransitionIndicators.jl",
    authors = "George Datseris",
    pages = DYNAMICALSYSTEMSBASE_PAGES,
    doctest = false,
    draft = false,
)

if CI
    deploydocs(
        repo = "github.com/JuliaDynamics/TransitionIndicators.jl.git",
        target = "build",
        push_preview = true
    )
end
