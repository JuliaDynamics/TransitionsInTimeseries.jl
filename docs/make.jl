cd(@__DIR__)

using TransitionsInTimeseries, Statistics, StatsBase

using Literate

# process examples and add then in a sidebar column
example_files = readdir(joinpath(@__DIR__, "src", "examples"))
example_pages = String[]
for file in example_files
    mkdownname = splitext(file)[1]*".md"
    Literate.markdown("src/examples/$(file)", "src/examples"; credit = false)
    push!(example_pages, "examples/$(mkdownname)")
end

pages = [
    "index.md",
    "tutorial.md",
    "api.md",
    "Examples" => example_pages,
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