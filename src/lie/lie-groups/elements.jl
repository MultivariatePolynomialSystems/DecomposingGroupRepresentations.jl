export GroupElem,
    matrix

struct GroupElem{T <: AbstractGroup, S <: AbstractMatrix} <: AbstractGroupElem
    group::T
    matrix::S
end

group(g::GroupElem) = g.group
matrix(g::GroupElem) = g.matrix

Base.rand(
    G::ScalingLieGroup{F}
) where F = GroupElem(G, prod([Diagonal(rand(F) .^ e) for e in eachrow(exponents(G))]))

Base.rand(
    G::LieGroup{<:Union{AbstractFloat, Complex{<:AbstractFloat}}}
) = GroupElem(G, exp(matrix(rand(algebra(G)))))
function Base.rand(G::LieGroup{F}) where F
    if name(G) == "SO(3)"
        println(F)
        return GroupElem(G, rand_rotation(F))
    end
end


struct DirectProductGroupElem{T <: AbstractDirectProductGroup, S <: GroupElem} <: AbstractGroupElem
    group::T
    elements::Vector{S}
end

group(g::DirectProductGroupElem) = g.group
elements(g::DirectProductGroupElem) = g.elements
Base.rand(G::DirectProductGroup) = DirectProductGroupElem(G, [rand(Gᵢ) for Gᵢ in groups(G)])

