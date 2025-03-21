```@meta
CurrentModule = DecomposingRepresentations
```

# Introduction
DecomposingRepresentations.jl is a Julia package that provides an API for decomposing vector space representations of reductive groups together with a concrete implementation of group actions on multivariate polynomials using [`DynamicPolynomials.jl`](https://github.com/JuliaAlgebra/DynamicPolynomials.jl).

## Quick start

```@repl
using DecomposingRepresentations

@polyvar x y z
vars = [x, y, z]

SO3 = LieGroup("SO", 3)
a = MatrixGroupAction(SO3, [vars])

V = FixedDegreePolynomials(vars, 2)
ρ = GroupRepresentation(a, V)

irr_decomp = IrreducibleDecomposition(ρ)
sum([dim(irr) for irr in irreducibles(irr_decomp)])

basis(space(irr_decomp[1])) # invariant
basis(space(irr_decomp[2])) # spherical harmonics
```
