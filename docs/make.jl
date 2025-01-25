using DecomposingRepresentations
using Documenter

DocMeta.setdocmeta!(DecomposingRepresentations, :DocTestSetup, :(using DecomposingRepresentations); recursive=true)

makedocs(;
    modules=[DecomposingRepresentations],
    authors="Viktor Korotynskiy <korotynskiy.viktor@gmail.com> and contributors",
    sitename="DecomposingRepresentations.jl",
    format=Documenter.HTML(;
        canonical="https://azoviktor.github.io/DecomposingRepresentations.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/azoviktor/DecomposingRepresentations.jl",
    devbranch="main",
)
