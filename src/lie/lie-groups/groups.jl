export LieGroup, DirectProductLieGroup, ×


struct ScalingLieGroup{F} <: AbstractReductiveLieGroup{F}
    name::String
    algebra::ScalingLieAlgebra{F}
end

ScalingLieGroup{F}(size::Int) where F = ScalingLieGroup("ℂˣ", ScalingLieAlgebra{F}(size))

name(G::ScalingLieGroup) = G.name
algebra(G::ScalingLieGroup) = G.algebra
exponents(G::ScalingLieGroup) = exponents(algebra(G))


struct LieGroup{F, T <: LieAlgebra{F}} <: AbstractReductiveLieGroup{F}
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


struct DirectProductLieGroup{F, T <: AbstractReductiveLieGroup{F}} <: AbstractDirectProductGroup{F}
    name::String
    groups::Vector{T}
end

DirectProductLieGroup(
    groups::Vector{AbstractReductiveLieGroup}
) = DirectProductLieGroup(join([name(G) for G in groups], " × "), groups)

name(G::DirectProductLieGroup) = G.name
groups(G::DirectProductLieGroup) = G.groups
algebra(G::DirectProductLieGroup) = SumLieAlgebra([algebra(Gᵢ) for Gᵢ in groups(G)])
×(
    G₁::AbstractReductiveLieGroup{F},
    G₂::AbstractReductiveLieGroup{F}
) where F = DirectProductLieGroup("$(name(G₁)) × $(name(G₂))", [G₁, G₂])

function Base.show(io::IO, G::DirectProductLieGroup{F}) where F
    println(io, "DirectProductLieGroup $(name(G))")
    println(io, " number type (or field): $(F)")
    print(io, " Lie algebra:")
end