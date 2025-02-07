export LieGroup, DirectProductLieGroup, ×


struct LieGroup{T <: AbstractReductiveLieAlgebra} <: AbstractReductiveLieGroup
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

function Base.show(io::IO, G::LieGroup)
    println(io, "LieGroup $(name(G))")
    println(io, " Lie algebra:")
    show(io, algebra(G); offset = 2)
end


struct DirectProductLieGroup <: AbstractReductiveLieGroup
    name::String
    groups::Vector{AbstractReductiveLieGroup}
end

DirectProductLieGroup(
    groups::Vector{AbstractReductiveLieGroup}
) = DirectProductLieGroup(join([name(G) for G in groups], " × "), groups)

name(G::DirectProductLieGroup) = G.name
groups(G::DirectProductLieGroup) = G.groups
algebra(G::DirectProductLieGroup) = SumLieAlgebra([algebra(Gᵢ) for Gᵢ in groups(G)])
×(G₁::AbstractReductiveLieGroup, G₂::AbstractReductiveLieGroup) = DirectProductLieGroup([G₁, G₂])

function Base.show(io::IO, G::DirectProductLieGroup)
    println(io, "DirectProductLieGroup $(name(G))")
    println(io, " Lie algebra:")
    show(io, algebra(G); offset = 2)
end