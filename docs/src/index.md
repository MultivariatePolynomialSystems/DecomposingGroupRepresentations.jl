```@meta
CurrentModule = DecomposingRepresentations
```

# Introduction
DecomposingRepresentations.jl is a Julia package that provides an API for decomposing vector space representations of reductive groups together with a concrete implementation of group actions on multivariate polynomials using [`DynamicPolynomials.jl`](https://github.com/JuliaAlgebra/DynamicPolynomials.jl).

## Quick start

```@repl
using DecomposingRepresentations

# ----- EXAMPLE 1 -----
@polyvar x y z
vars = [x, y, z]

SO3 = LieGroup("SO", 3)
a = MatrixGroupAction(SO3, [vars])

V = FixedDegreePolynomials(vars, 2)
ρ = GroupRepresentation(a, V)

irr_decomp = IrreducibleDecomposition(ρ)
sum([dim(irr) for irr in irreducibles(irr_decomp)])

basis(space(irr_decomp[1]))
basis(space(irr_decomp[2]))

# ----- EXAMPLE 2 -----
@polyvar x[1:3] y[1:3]

SO3 = LieGroup("SO", 3)
a = MatrixGroupAction(SO3, [x, y])

V = FixedDegreePolynomials(vcat(x, y), 2)
ρ = GroupRepresentation(a, V)

irr_decomp = IrreducibleDecomposition(ρ)
sum([dim(irr) for irr in irreducibles(irr_decomp)])

iso_decomp = IsotypicDecomposition(ρ)
inv = invariant_component(iso_decomp)
basis(space(inv))

# ----- EXAMPLE 3 -----
@polyvar x[1:3] y[1:3]

SO3 = LieGroup("SO", 3)
a = MatrixGroupAction(SO3, [x, y])

V = FixedMultidegreePolynomials([x, y], [1, 1])
ρ = GroupRepresentation(a, V)

irr_decomp = IrreducibleDecomposition(ρ)
sum([dim(irr) for irr in irreducibles(irr_decomp)])

inv = invariant_component(IsotypicDecomposition(ρ))
basis(space(inv))
```
