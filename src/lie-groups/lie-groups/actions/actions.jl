export MatrixGroupAction,
    action_vectors,
    ScalingLieGroupAction,
    action_vector,
    DirectProductGroupAction,
    ×,
    are_commutative,
    act

struct MatrixGroupAction{T<:GroupType, F, S<:AbstractGroup{T,F}, V<:Variable} <: AbstractGroupAction{T, F}
    group::S
    vars::Vector{Vector{V}}
end

MatrixGroupAction(
    G::S,
    vectors::AbstractVector{<:AbstractVector{V}}
) where {T<:GroupType, F, S <: AbstractGroup{T, F}, V <: Variable} = MatrixGroupAction{T, F, S, V}(G, collect(collect(v) for v in vectors))

group(a::MatrixGroupAction) = a.group
action_vectors(a::MatrixGroupAction) = a.vars
DynamicPolynomials.variables(a::MatrixGroupAction) = vcat(action_vectors(a)...)
space(a::MatrixGroupAction{T, F}) where {T,F} = VectorSpace(F, variables(a))

function Base.show(io::IO, a::MatrixGroupAction)
    println(io, "MatrixGroupAction of $(name(group(a)))")
    print(io, " vectors under action: ")
    print(io, "[", join([join(map(repr, v), ", ") for v in action_vectors(a)], "], ["), "]")
end


# TODO: remove and leave MatrixGroupAction only?
struct ScalingLieGroupAction{F, T <: ScalingLieGroup{F}, V <: Variable} <: AbstractGroupAction{Lie, F}
    group::T
    vars::Vector{V}
end

group(a::ScalingLieGroupAction) = a.group
action_vector(a::ScalingLieGroupAction) = a.vars
DynamicPolynomials.variables(a::ScalingLieGroupAction) = action_vector(a)
space(a::ScalingLieGroupAction{F}) where F = VectorSpace(F, action_vector(a))

ScalingLieGroupAction(v::Vector{<:Variable}) = ScalingLieGroupAction(ScalingLieGroup{ComplexF64}(length(v)), v)

function ScalingLieGroupAction(V::AbstractMatrix{<:Variable})
    exps = spzeros(Int, size(V, 2), length(V))
    for j in 1:size(V, 2)
        exps[j, ((j-1)*size(V,1)+1):(j*size(V,1))] = ones(Int, size(V, 1))
    end
    F = ComplexF64
    alg = ScalingLieAlgebra{F}("(ℂ)$(superscript(size(V,2)))", exps)
    G = ScalingLieGroup{F}("(ℂˣ)$(superscript(size(V,2)))", alg)
    return ScalingLieGroupAction(G, V[:])
end

MatrixGroupAction(a::ScalingLieGroupAction) = MatrixGroupAction(group(a), [action_vector(a)])

function show_action(io::IO, a::ScalingLieGroupAction; offset::Int=0)
    U = exponents(group(a))
    if size(U, 1) == 1
        @polyvar λ
        λ = [λ]
    else
        @polyvar λ[1:size(U, 1)]
    end
    action = []
    vars = action_vector(a)
    for j in axes(U, 1)
        nzind, nzval = findnz(U[j, :])
        exprs = (λ[j].^nzval).*vars[nzind]
        push!(action, collect(zip(vars[nzind], exprs)))
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
    println(io, " vector under action: ", "[", join(map(repr, action_vector(a)), ", "), "]")
    print(io, " action:")
    show_action(io, a; offset = 2)
end

struct DirectProductGroupAction{T<:GroupType, F, S<:AbstractDirectProductGroup{T, F}} <: AbstractGroupAction{T, F}
    group::S
    lie_actions::Vector{AbstractGroupAction{Lie, F}}
    finite_actions::Vector{AbstractGroupAction{Finite, F}}
end

DirectProductGroupAction(
    G::T,
    lie_actions::AbstractVector{S₁}
) where {
    F,
    T <: AbstractDirectProductGroup{Lie, F},
    S₁ <: AbstractGroupAction{Lie, F}
} = DirectProductGroupAction{Lie, F, T}(G, lie_actions, [])

group(a::DirectProductGroupAction) = a.group
lie_actions(a::DirectProductGroupAction) = a.lie_actions
finite_actions(a::DirectProductGroupAction) = a.finite_actions
actions(a::DirectProductGroupAction) = vcat(lie_actions(a), finite_actions(a))
nactions(a::DirectProductGroupAction) = length(lie_actions(a)) + length(finite_actions(a))
space(a::DirectProductGroupAction) = +([space(a) for a in actions(a)]...)

function Base.show(io::IO, a::DirectProductGroupAction)
    println(io, "DirectProductGroupAction of $(name(group(a)))")
    print(io, " lie actions:")
    # for action in lie_actions(a)
    #     println(io)
    #     show_action(io, action; offset=2, show_name=true)
    # end
end

function are_commutative(
    a₁::AbstractGroupAction{T₁, F},
    a₂::AbstractGroupAction{T₂, F}
) where {T₁<:GroupType, T₂<:GroupType, F}
    V = space(a₁) + space(a₂)
    v = rand(V)
    g₁, g₂ = rand(group(a₁)), rand(group(a₂))
    v₁ = act(g₂, a₂, act(g₁, a₁, v))
    v₂ = act(g₁, a₁, act(g₂, a₂, v))
    return v₁ ≈ v₂
end

function ×(
    a₁::AbstractGroupAction{Lie},
    a₂::AbstractGroupAction{Lie}
)
    @assert are_commutative(a₁, a₂)
    return DirectProductGroupAction(group(a₁) × group(a₂), [a₁, a₂])
end

function ×(
    a₁::DirectProductGroupAction{Lie},
    a₂::AbstractGroupAction{Lie}
)
    @assert are_commutative(a₁, a₂)
    return DirectProductGroupAction(group(a₁) × group(a₂), vcat(lie_actions(a₁), [a₂]))
end