using Backboner
using Documenter

DocMeta.setdocmeta!(Backboner, :DocTestSetup, :(using Backboner, Backboner.Protein); recursive=true)

makedocs(;
    modules = [Backboner, Backboner.Protein],
    sitename = "Backboner.jl",
    format = Documenter.HTML(prettyurls = (get(ENV, "CI", "false") == "true")),
    doctest = false,
    pages = [
        "Overview" => "index.md",
        "Backboner API" => "backboner.md",
        "Protein API" => "protein.md",
    ],
    authors = "Anton Oresten",
    checkdocs = :all,
)

deploydocs(;
    repo = "github.com/MurrellGroup/Backboner.jl.git",
    push_preview = true,
    devbranch = "main",
    deps = nothing,
    make = nothing,
)
