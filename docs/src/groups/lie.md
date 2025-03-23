# Lie groups

The computations with Lie groups are done by passing to their associated Lie algebras. We distinguish between Lie groups that act by scalings (represented by [`ScalingLieGroup`](@ref)) and basic reductive matrix Lie groups, like ``\mathrm{SO}(n)`` (represented by [`LieGroup`](@ref)).

## Lie groups

```@docs
ScalingLieGroup
LieGroup
```

## Lie algebras

```@docs
AbstractLieAlgebra
algebra(::AbstractGroup{Lie})
name(::AbstractLieAlgebra)
basis(::AbstractLieAlgebra)
chevalley_basis(::AbstractLieAlgebra)
dim(::AbstractLieAlgebra)
rank(::AbstractLieAlgebra)
```