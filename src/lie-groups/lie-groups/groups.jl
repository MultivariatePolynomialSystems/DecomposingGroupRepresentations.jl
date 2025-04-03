export LieGroup,
    ScalingLieGroup,
    DirectProductGroup,
    ×


"""
    ScalingLieGroup{F} <: AbstractGroup{Lie, F}

Represents a scaling Lie group. The group elements are diagonal matrices.

# Constructors
```julia
ScalingLieGroup{F}(size::Int) where F
ScalingLieGroup{F}(exps::Matrix{Int}) where F
```

# Examples
```jldoctest
julia> ScalingLieGroup{ComplexF64}(5)
ScalingLieGroup ℂˣ
 number type (or field): ComplexF64
 Lie algebra properties:
  dimension: 1
  basis (diagonal matrices):
   [1, 1, 1, 1, 1]
```
```jldoctest
julia> exps = [1 1 1 0 0 0; 0 0 0 1 -2 0];

julia> ScalingLieGroup{ComplexF64}(exps)
ScalingLieGroup (ℂˣ)²
 number type (or field): ComplexF64
 Lie algebra properties:
  dimension: 2
  basis (diagonal matrices):
   [1, 1, 1, 0, 0, 0]
   [0, 0, 0, 1, -2, 0]
```
"""
struct ScalingLieGroup{F} <: AbstractGroup{Lie, F}
    name::String
    algebra::ScalingLieAlgebra{F}
end

ScalingLieGroup{F}(size::Int) where F = ScalingLieGroup("ℂˣ", ScalingLieAlgebra{F}(size))
ScalingLieGroup(size::Int) = ScalingLieGroup{ComplexF64}(size)
function ScalingLieGroup{F}(exps::Matrix{Int}) where F
    name = size(exps, 1) == 1 ? "ℂˣ" : "(ℂˣ)$(superscript(size(exps, 1)))"
    ScalingLieGroup(name, ScalingLieAlgebra{F}(exps))
end
ScalingLieGroup(exps::Matrix{Int}) = ScalingLieGroup{ComplexF64}(exps)

name(G::ScalingLieGroup) = G.name
algebra(G::ScalingLieGroup) = G.algebra
exponents(G::ScalingLieGroup) = exponents(algebra(G))
weight(G::ScalingLieGroup, i::Integer) = Weight(Vector(exponents(G)[:, i]))

function Base.show(io::IO, G::ScalingLieGroup{F}) where F
    println(io, "ScalingLieGroup $(name(G))")
    println(io, " number type (or field): $(F)")
    println(io, " Lie algebra properties:")
    println(io, "  dimension: $(dim(algebra(G)))")
    println(io, "  basis (diagonal matrices):")
    show_basis(io, algebra(G), offset=3)
end


"""
    LieGroup{F} <: AbstractGroup{Lie, F}

Describes a matrix Lie group.

# Constructors
```julia
LieGroup(type::String, size::Int)
```

Supported Lie group types: SO (special orthogonal).

# Examples
```jldoctest
julia> SO3 = LieGroup("SO", 3)
LieGroup SO(3)
 number type (or field): ComplexF64
 weight type: Int64
 Lie algebra properties:
  dimension: 3
  rank (dimension of Cartan subalgebra): 1
```
"""
struct LieGroup{F, T <: LieAlgebra{F}} <: AbstractGroup{Lie, F}
    name::String
    algebra::T
end

name(G::LieGroup) = G.name
algebra(G::LieGroup) = G.algebra

SO3() = LieGroup("SO(3)", so3())

# TODO
function LieGroup(name::String, size::Int)
    if name == "SO" && size == 3
        return SO3()
    else
        error("Not implemented")
    end
end

function Base.show(io::IO, G::LieGroup{F}) where F
    println(io, "LieGroup $(name(G))")
    println(io, " number type (or field): $(F)")
    println(io, " weight type: $(weight_inner_type(algebra(G)))")
    println(io, " Lie algebra properties:")
    println(io, "  dimension: $(dim(algebra(G)))")
    print(io, "  rank (dimension of Cartan subalgebra): $(rank(algebra(G)))")
end

"""
    DirectProductGroup{T<:GroupType, F} <: AbstractDirectProductGroup{T, F}

Represents a direct product of groups.
# Examples
```jldoctest
julia> SO3 = LieGroup("SO", 3);

julia> T = ScalingLieGroup{ComplexF64}([1 2 3 4; -1 -2 -3 -4]);

julia> SO3 × SO3 × T
DirectProductGroup SO(3) × SO(3) × (ℂˣ)²
 number type (or field): ComplexF64
 3 factors: SO(3), SO(3), (ℂˣ)²
 Lie algebra:
  SumLieAlgebra 𝖘𝖔(3) ⊕ 𝖘𝖔(3) ⊕ ℂ²
  dimension: 8
  rank (dimension of Cartan subalgebra): 4
```
"""
struct DirectProductGroup{T<:GroupType, F, S<:AbstractGroup{T, F}} <: AbstractDirectProductGroup{T, F}
    name::String
    groups::Vector{S}
end

DirectProductGroup(
    groups::Vector{<:AbstractGroup{Lie}}
) = DirectProductGroup(join([name(G) for G in groups], " × "), groups)

name(G::DirectProductGroup) = G.name
groups(G::DirectProductGroup) = G.groups
ngroups(G::DirectProductGroup) = length(groups(G))
algebra(G::DirectProductGroup) = SumLieAlgebra([algebra(Gᵢ) for Gᵢ in groups(G)])
×(
    G₁::AbstractGroup{Lie, F},
    G₂::AbstractGroup{Lie, F}
) where F = DirectProductGroup("$(name(G₁)) × $(name(G₂))", [G₁, G₂])
×(
    G₁::DirectProductGroup{Lie, F},
    G₂::AbstractGroup{Lie, F}
) where F = DirectProductGroup("$(name(G₁)) × $(name(G₂))", vcat(groups(G₁), [G₂]))

function Base.show(io::IO, G::DirectProductGroup{T, F}) where {T, F}
    println(io, "DirectProductGroup $(name(G))")
    println(io, " number type (or field): $(F)")
    println(io, " $(ngroups(G)) factors: ", join([name(Gᵢ) for Gᵢ in groups(G)], ", "))
    println(io, " Lie algebra:")
    println(io, "  SumLieAlgebra $(name(algebra(G)))")
    println(io, "  dimension: $(dim(algebra(G)))")
    print(io, "  rank (dimension of Cartan subalgebra): $(rank(algebra(G)))")
end

# TODO: remove and use DirectProductGroup instead with a Mixed group type
struct DirectProductMixedGroup{
    F, T <: AbstractGroup{Lie, F}, S <: AbstractGroup{Finite, F}
} <: AbstractDirectProductGroup{Mixed, F}
    name::String
    lie_part::T
    finite_part::S
end

zero_weight(G::AbstractGroup{Lie}) = zero_weight(algebra(G))