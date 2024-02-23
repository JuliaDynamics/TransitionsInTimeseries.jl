cd(@__DIR__)

using TransitionsInTimeseries, StatsBase
using Literate

# process examples and add then in a sidebar column
jl_examples = ["logistic.jl", "ks_paleojump.jl", "do-events.jl"]
example_pages = String[]
for file in jl_examples
    mkdownname = splitext(file)[1]*".md"
    Literate.markdown("src/examples/$(file)", "src/examples"; credit = false)
    push!(example_pages, "examples/$(mkdownname)")
end
Literate.markdown("src/tutorial.jl", "src"; credit = false)

pages = [
    "index.md",
    "tutorial.md",
    "Examples" => example_pages,
    "api.md",
    "refs.md",
    "devdocs.md",
]

import Downloads
Downloads.download(
    "https://raw.githubusercontent.com/JuliaDynamics/doctheme/master/build_docs_with_style.jl",
    joinpath(@__DIR__, "build_docs_with_style.jl")
)
include("build_docs_with_style.jl")

using DocumenterCitations
bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:authoryear)

build_docs_with_style(pages, TransitionsInTimeseries, StatsBase;
    authors = "Jan Swierczek-Jereczek <jan.jereczek@gmail.com>, "*
    "George Datseris <datseris.george@gmail.com>", bib
)
