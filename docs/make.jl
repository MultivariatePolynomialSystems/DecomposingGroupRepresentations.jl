using DecomposingGroupRepresentations
using Documenter

DocMeta.setdocmeta!(DecomposingGroupRepresentations, :DocTestSetup, :(using DecomposingGroupRepresentations); recursive=true)

makedocs(;
    modules=[DecomposingGroupRepresentations],
    authors="Viktor Korotynskiy <korotynskiy.viktor@gmail.com> and contributors",
    repo="https://github.com/MultivariatePolynomialSystems/DecomposingGroupRepresentations.jl",
    sitename="DecomposingGroupRepresentations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://multivariatepolynomialsystems.github.io/DecomposingGroupRepresentations.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "Reductive groups" => [
            "Types of groups" => "groups/types.md",
            "Finite groups" => "groups/finite.md",
            "Lie groups" => "groups/lie.md",
            "Direct products" => "groups/products.md",
        ],
        "Representations of reductive groups" => [
            "Vector spaces" => "representations/spaces.md",
            "Actions" => "representations/actions.md",
            "Representations" => "representations/reprs.md",
        ],
        "Decomposing representations" => [
            "Irreducible decomposition" => "decompose/irreducible.md",
            "Isotypic decomposition" => "decompose/isotypic.md",
        ],
    ],
)

deploydocs(;
    repo="github.com/MultivariatePolynomialSystems/DecomposingGroupRepresentations.jl.git",
    devbranch="main",
    versions = ["stable", "dev"]
)
