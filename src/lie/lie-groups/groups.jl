export ReductiveLieGroup


struct ReductiveLieGroup{T <: AbstractReductiveLieAlgebra} <: AbstractReductiveLieGroup
    name::String
    algebra::T
end

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
    println(io, "ReductiveLieGroup $(G.name)")
    print(io, " Lie algebra: ")
end
