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
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/MultivariatePolynomialSystems/MultivariateInterpolation.jl.git",
    devbranch="main",
)
