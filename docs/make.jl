using Backboner
using Documenter

DocMeta.setdocmeta!(
    Backboner,
    :DocTestSetup,
    :(using Backboner);
    recursive=true,
)

makedocs(;
    modules=[Backboner],
    authors="Anton Oresten <antonoresten@gmail.com> and contributors",
    sitename="Backboner.jl",
    format=Documenter.HTML(;
        canonical="https://MurrellGroup.github.io/Backboner.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "API.md",
    ],
    doctest=true,
)

deploydocs(;
    repo="github.com/MurrellGroup/Backboner.jl",
    devbranch="main",
)