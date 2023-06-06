cd(@__DIR__)

using TransitionsInTimeseries, Statistics, StatsBase

pages = [
    "index.md",
    "tutorial.md",
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