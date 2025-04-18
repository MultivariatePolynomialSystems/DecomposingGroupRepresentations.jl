# Lie groups
Lie groups in this Julia package are represented by the abstract type
```julia
AbstractGroup{Lie, F}
```
where `F` defines the field (or number type) over which the group is defined.

The computations with Lie groups are done by passing to their associated Lie algebras. We distinguish between Lie groups that act by scalings (represented by [`ScalingLieGroup`](@ref)) and basic reductive matrix Lie groups, like ``\mathrm{SO}(n)`` (represented by [`LieGroup`](@ref)).

## Lie groups

```@docs
ScalingLieGroup
LieGroup
```

!!! note
    Supported Lie group types and sizes: SO(3).

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