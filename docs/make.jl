cd(@__DIR__)

using TransitionsInTimeseries, Statistics, StatsBase, Literate

Literate.markdown("src/examples/tutorial.jl", "src/examples"; credit = false)
# Literate.markdown("src/examples/permutation_entropy.jl", "src"; credit = false)
Literate.markdown("src/examples/do-events.jl", "src/examples"; credit = false)

pages = [
    "index.md",
    "examples/tutorial.md",
    # "examples/permutation_entropy.md",
    "examples/do-events.md",
    "api.md",
]

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

build_docs_with_style(pages, TransitionsInTimeseries, Statistics, StatsBase;
    authors = "Jan Swierczek-Jereczek <jan.jereczek@gmail.com>, George Datseris <datseris.george@gmail.com>"
)