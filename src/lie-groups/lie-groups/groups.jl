export LieGroup,
    ScalingLieGroup,
    DirectProductGroup,
    ×


struct ScalingLieGroup{F} <: AbstractGroup{Lie, F}
    name::String
    algebra::ScalingLieAlgebra{F}
end

ScalingLieGroup{F}(size::Int) where F = ScalingLieGroup("ℂˣ", ScalingLieAlgebra{F}(size))

name(G::ScalingLieGroup) = G.name
algebra(G::ScalingLieGroup) = G.algebra
exponents(G::ScalingLieGroup) = exponents(algebra(G))
weight(G::ScalingLieGroup, i::Integer) = Weight(Vector(exponents(G)[:, i]))


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
    println(io, " Lie algebra:")
    show(io, algebra(G); offset = 2)
end


struct DirectProductGroup{T<:GroupType, F, S<:AbstractGroup{T, F}} <: AbstractDirectProductGroup{T, F}
    name::String
    groups::Vector{S}
end

DirectProductGroup(
    groups::Vector{<:AbstractGroup{Lie}}
) = DirectProductGroup(join([name(G) for G in groups], " × "), groups)

name(G::DirectProductGroup) = G.name
groups(G::DirectProductGroup) = G.groups
algebra(G::DirectProductGroup) = SumLieAlgebra([algebra(Gᵢ) for Gᵢ in groups(G)])
×(
    G₁::AbstractGroup{Lie, F},
    G₂::AbstractGroup{Lie, F}
) where F = DirectProductGroup("$(name(G₁)) × $(name(G₂))", [G₁, G₂])
×(
    G₁::DirectProductGroup{Lie, F},
    G₂::AbstractGroup{Lie, F}
) where F = DirectProductGroup("$(name(G₁)) × $(name(G₂))", vcat(groups(G₁), [G₂]))

function Base.show(io::IO, G::DirectProductGroup{F}) where F
    println(io, "DirectProductLieGroup $(name(G))")
    println(io, " number type (or field): $(F)")
    print(io, " Lie algebra:")
end


struct DirectProductMixedGroup{
    F, T <: AbstractGroup{Lie, F}, S <: AbstractGroup{Finite, F}
} <: AbstractDirectProductGroup{Mixed, F}
    name::String
    lie_part::T
    finite_part::S
end

zero_weight(G::AbstractGroup{Lie}) = zero_weight(algebra(G))