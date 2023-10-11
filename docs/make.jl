cd(@__DIR__)

using TransitionsInTimeseries, Statistics, StatsBase
using Literate

# process examples and add then in a sidebar column
example_files = readdir(joinpath(@__DIR__, "src", "examples"))
jl_indices = [f[end-2:end] == ".jl" for f in example_files]
jl_examples = example_files[jl_indices]

example_pages = String[]
for file in jl_examples
    mkdownname = splitext(file)[1]*".md"
    Literate.markdown("src/examples/$(file)", "src/examples"; credit = false)
    push!(example_pages, "examples/$(mkdownname)")
end

# Sort pages with increasing complexity rather than alphabetically
permute!(example_pages, [2, 1])

pages = [
    "index.md",
    "tutorial.md",
    "Examples" => example_pages,
    "api.md",
    "refs.md",
]

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

using DocumenterCitations
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)

build_docs_with_style(pages, TransitionsInTimeseries, Statistics, StatsBase;
    authors = "Jan Swierczek-Jereczek <jan.jereczek@gmail.com>, "*
    "George Datseris <datseris.george@gmail.com>", bib)