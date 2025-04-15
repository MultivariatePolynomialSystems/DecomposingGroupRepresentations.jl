```@meta
CurrentModule = DecomposingGroupRepresentations
```

# Introduction
DecomposingGroupRepresentations.jl is a Julia package that provides an API for decomposing representations of reductive groups acting on multivariate polynomials using [`DynamicPolynomials.jl`](https://github.com/JuliaAlgebra/DynamicPolynomials.jl).

## Quick start

```@repl
using DecomposingGroupRepresentations

@polyvar x y z
vars = [x, y, z]

SO3 = LieGroup("SO", 3)
a = MatrixGroupAction(SO3, [vars])

V = FixedDegreePolynomials(vars, 2)
ρ = GroupRepresentation(a, V)

irrs = IrreducibleDecomposition(ρ)
[highest_weight(irr) for irr in irreducibles(irrs)]
[vector(hw_vector(irr)) for irr in irreducibles(irrs)]

iso = IsotypicDecomposition(ρ)
basis(iso[Weight([0])]) # invariant
H₂ = basis(iso[Weight([2])]) # spherical harmonics
rref(H₂)
```
