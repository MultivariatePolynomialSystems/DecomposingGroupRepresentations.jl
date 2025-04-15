# DecomposingRepresentations.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://multivariatepolynomialsystems.github.io/DecomposingRepresentations.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-purple.svg)](https://multivariatepolynomialsystems.github.io/DecomposingRepresentations.jl/dev/)
<!-- [![Build Status](https://github.com/multivariatepolynomialsystems/DecomposingRepresentations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/multivariatepolynomialsystems/DecomposingRepresentations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/multivariatepolynomialsystems/DecomposingRepresentations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/azoviktor/DecomposingRepresentations.jl) -->

DecomposingRepresentations.jl is a Julia package for decomposing representations of reductive groups acting on multivariate polynomials.

## Installation

Enter the Pkg REPL by pressing `]` from the Julia REPL and then type
```julia
add DecomposingRepresentations
```
To get back to the Julia REPL, press backspace.

## Usage
```julia
using DecomposingRepresentations

@polyvar x y z
vars = [x, y, z]
SO3 = LieGroup("SO", 3)
a = MatrixGroupAction(SO3, [vars])
V = FixedDegreePolynomials(vars, 2)
ρ = GroupRepresentation(a, V)

irrs = irreducibles(ρ)
```
```
IrreducibleDecomposition of SO(3)-action on 6-dimensional vector space
 number of irreducibles: 2
 dimensions of irreducibles: 1, 5
```