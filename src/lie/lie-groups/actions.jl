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
exponents(a::ScalingLieGroupAction) = exponents(algebra(group(a)))

ScalingLieGroupAction(v::Vector) = ScalingLieGroupAction(ReductiveLieGroup("ℂˣ", ScalingLieAlgebra(length(v))), v)

function show_action(io::IO, a::ScalingLieGroupAction{<:Variable}; offset::Int=0)
    U = exponents(a)
    if size(U, 1) == 1
        @polyvar λ
        λ = [λ]
    else
        @polyvar λ[1:size(U, 1)]
    end
    action = Vector{Vector{Tuple{Variable, Monomial}}}([[] for _ in axes(U, 1)])
    vars = action_vector(a)
    for j in axes(U, 1)
        nzind, nzval = findnz(U[j, :])
        exprs = (λ[j].^nzval).*vars[nzind]
        action[j] = collect(zip(vars[nzind], exprs))
    end
    for free_action in action
        print(io, "\n", " "^offset)
        for (j, (var, expr)) in enumerate(free_action)
            print(io, repr(var), " ↦ ", repr(expr))
            j < length(free_action) && print(io, ", ")
        end
    end
end

function Base.show(io::IO, a::ScalingLieGroupAction)
    println(io, "ScalingLieGroupAction of $(name(group(a)))")
    print(io, " action:")
    show_action(io, a; offset = 2)
end