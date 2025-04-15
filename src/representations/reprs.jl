export GroupRepresentation,
    irreducibles,
    common_nullspace_as_weight_vectors, ws_nullspace,
    isotypics


"""
    GroupRepresentation <: AbstractGroupRepresentation

Represents a group representation. Consists of a group action and a vector space.

# Constructors
```julia
GroupRepresentation(a::AbstractGroupAction, V::AbstractSpace)
```

# Examples
```jldoctest
julia> @polyvar x[1:3];

julia> SO3 = LieGroup("SO", 3);

julia> a = MatrixGroupAction(SO3, [x]);

julia> V = FixedDegreePolynomials(x, 2);

julia> ρ = GroupRepresentation(a, V)
GroupRepresentation of SO(3) on 6-dimensional vector space
 Lie group: SO(3)

julia> dim(ρ)
6

julia> group(ρ)
LieGroup SO(3)
 number type (or field): ComplexF64
 weight type: Int64
 Lie algebra properties:
  dimension: 3
  rank (dimension of Cartan subalgebra): 1
```
"""
struct GroupRepresentation{A<:AbstractGroupAction{Lie}, T<:AbstractSpace} <: AbstractGroupRepresentation{Lie, T}
    action::A
    V::T
end

action(ρ::GroupRepresentation) = ρ.action
space(ρ::GroupRepresentation) = ρ.V

# called by Shift+Enter
function Base.show(io::IO, ::MIME"text/plain", ρ::GroupRepresentation)
    println(
        io,
        "GroupRepresentation of $(name(group(ρ))) ",
        "on $(dim(ρ))-dimensional vector space"
    )
    print(io, " Lie group: ", name(group(ρ)))
end

# called by print and inside vectors/matrices
function Base.show(io::IO, ρ::GroupRepresentation)
    print(
        io,
        "GroupRepresentation of $(name(group(ρ))) ",
        "on $(dim(ρ))-dimensional vector space"
    )
end

summands(
    ρ::GroupRepresentation{<:AbstractGroupAction, <:AbstractDirectSum}
) = [GroupRepresentation(action(ρ), V) for V in summands(space(ρ))]

GroupRepresentation(a::AbstractGroupAction{Lie}, V::HighestWeightModule) = IrreducibleRepresentation(a, V)