# Groups

```@docs
AbstractGroup
```

## Group types

```@docs
GroupType
Finite
Lie
Mixed
```

## Finite groups

## Lie groups

```@docs
algebra(::AbstractGroup{Lie})
```

### Lie algebras

```@docs
AbstractLieAlgebra
name(::AbstractLieAlgebra)
basis(::AbstractLieAlgebra)
chevalley_basis(::AbstractLieAlgebra)
dim(::AbstractLieAlgebra)
rank(::AbstractLieAlgebra)
```

## Direct product of groups

```@docs
AbstractDirectProductGroup
```