export ReductiveLieGroup


struct ReductiveLieGroup{T <: AbstractReductiveLieAlgebra} <: AbstractReductiveLieGroup
    name::String
    algebra::T
end

name(G::ReductiveLieGroup) = G.name
algebra(G::ReductiveLieGroup) = G.algebra

SO3() = ReductiveLieGroup("SO(3)", so3())

# TODO
function ReductiveLieGroup(name::String, size::Int)
    if name == "SO" && size == 3
        return SO3()
    else
        error("Not implemented")
    end
end

function Base.show(io::IO, G::ReductiveLieGroup)
    println(io, "ReductiveLieGroup $(name(G))")
    println(io, " Lie algebra: ")
    show(io, algebra(G); offset = 2)
end
