```@meta
CurrentModule = DecomposingRepresentations
```

# Introduction
DecomposingRepresentations.jl is a Julia package that provides an API for decomposing vector space representations of reductive groups together with a concrete implementation of group actions on multivariate polynomials using [`DynamicPolynomials.jl`](https://github.com/JuliaAlgebra/DynamicPolynomials.jl).

## Quick start

```@repl
using DecomposingRepresentations
@polyvar x[1:3] y[1:3]
SO3 = LieGroup("SO", 3)
a = MatrixGroupAction(SO3, [x, y])
V = SymmetricPowersSpace(vcat(x, y), 0:2)
ρ = GroupRepresentation(a, V)
irrs = irreducibles(ρ)
sum([dim(irr) for irr in irrs])
iso = isotypic_components(ρ)
iso[2]
basis(space(iso[2]))
```
