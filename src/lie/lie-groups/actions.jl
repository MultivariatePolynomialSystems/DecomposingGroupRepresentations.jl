export MatrixGroupAction,
    action_vectors,
    ScalingLieGroupAction

struct MatrixGroupAction{T <: AbstractReductiveLieGroup, S} <: AbstractGroupAction
    group::T
    vectors::Vector{Vector{S}}
end

MatrixGroupAction(
    G::AbstractReductiveLieGroup,
    vectors
) = MatrixGroupAction(G, Vector{Vector{eltype(first(vectors))}}(vectors))

group(a::MatrixGroupAction) = a.group
action_vectors(a::MatrixGroupAction) = a.vectors

function Base.show(io::IO, a::MatrixGroupAction)
    println(io, "MatrixGroupAction of $(name(group(a)))")
    print(io, " vectors under action: ")
    print(io, "[", join([join(map(repr, v), ", ") for v in action_vectors(a)], "], ["), "]")
end


struct ScalingLieGroupAction{T} <: AbstractGroupAction
    group::ReductiveLieGroup{ScalingLieAlgebra}
    vector::Vector{T}
end

group(a::ScalingLieGroupAction) = a.group
action_vector(a::ScalingLieGroupAction) = a.vector
action_vectors(a::ScalingLieGroupAction) = [a.vector]

ScalingLieGroupAction(v::Vector) = ScalingLieGroupAction(ReductiveLieGroup("ℂˣ", ScalingLieAlgebra(length(v))), v)

function Base.show(io::IO, a::ScalingLieGroupAction)
    println(io, "ScalingLieGroupAction of $(name(group(a)))")
    print(io, " action:")
end