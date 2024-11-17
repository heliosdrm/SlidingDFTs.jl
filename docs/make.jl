using SlidingDFTs
using Documenter

DocMeta.setdocmeta!(SlidingDFTs, :DocTestSetup, :(using SlidingDFTs); recursive=true)

makedocs(;
    modules=[SlidingDFTs],
    authors="Helios De Rosario",
    sitename="SlidingDFTs.jl",
    format=Documenter.HTML(;
        canonical="https://heliosdrm.github.io/SlidingDFTs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Implementing SlidingDFTs" => "implement.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/heliosdrm/SlidingDFTs.jl",
    devbranch="main",
)
