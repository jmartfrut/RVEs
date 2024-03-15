using MimosaRVEs
using Documenter

DocMeta.setdocmeta!(MimosaRVEs, :DocTestSetup, :(using MimosaRVEs); recursive=true)

makedocs(;
    modules=[MimosaRVEs],
    authors="MultiSimo_Group",
    sitename="MimosaRVEs.jl",
    format=Documenter.HTML(;
        canonical="https://jmartfrut.github.io/MimosaRVEs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jmartfrut/MimosaRVEs.jl",
    devbranch="main",
)
