using DecomposingRepresentations
using Documenter

DocMeta.setdocmeta!(DecomposingRepresentations, :DocTestSetup, :(using DecomposingRepresentations); recursive=true)

makedocs(;
    modules=[DecomposingRepresentations],
    authors="Viktor Korotynskiy <korotynskiy.viktor@gmail.com> and contributors",
    repo="https://github.com/MultivariatePolynomialSystems/DecomposingRepresentations.jl",
    sitename="DecomposingRepresentations.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://multivariatepolynomialsystems.github.io/DecomposingRepresentations.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "Reductive groups" => [
            "Finite groups" => "groups/finite.md",
            "Lie groups" => "groups/lie.md",
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
    repo="github.com/MultivariatePolynomialSystems/DecomposingRepresentations.jl.git",
    devbranch="main",
)
