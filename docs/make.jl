cd(@__DIR__)

using TransitionIndicators

pages = [
    "index.md",
]

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

build_docs_with_style(pages, TransitionIndicators;
    authors = ["Jan Swierczek-Jereczek <jan.jereczek@gmail.com>", "George Datseris <datseris.george@gmail.com>"]
)