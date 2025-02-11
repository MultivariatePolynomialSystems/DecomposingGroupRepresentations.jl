export MatrixGroupAction,
    group,
    action_vectors,
    ScalingLieGroupAction,
    action_vector,
    ⊠,
    are_commutative

struct MatrixGroupAction{T <: AbstractGroup, S <: Variable} <: AbstractAction{T}
    group::T
    vars::Vector{Vector{S}}
end

MatrixGroupAction(
    G::T,
    vectors::AbstractVector{<:AbstractVector{S}}
) where {T <: AbstractGroup, S <: Variable} = MatrixGroupAction{T, S}(G, collect(collect(v) for v in vectors))

group(a::MatrixGroupAction) = a.group
action_vectors(a::MatrixGroupAction) = a.vars
space(a::MatrixGroupAction{<:AbstractGroup{T, F}}) where {T, F} = VariableSpace{F}(vcat(action_vectors(a)...))

function Base.show(io::IO, a::MatrixGroupAction)
    println(io, "MatrixGroupAction of $(name(group(a)))")
    print(io, " vectors under action: ")
    print(io, "[", join([join(map(repr, v), ", ") for v in action_vectors(a)], "], ["), "]")
end


struct ScalingLieGroupAction{T <: ScalingLieGroup, S <: Variable} <: AbstractAction{T}
    group::T
    vars::Vector{S}
end

group(a::ScalingLieGroupAction) = a.group
action_vector(a::ScalingLieGroupAction) = a.vars
space(a::ScalingLieGroupAction{ScalingLieGroup{F}}) where F = VariableSpace{F}(action_vector(a))

ScalingLieGroupAction(v::Vector{<:Variable}) = ScalingLieGroupAction(ScalingLieGroup{ComplexF64}(length(v)), v)
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


struct DirectProductGroupAction{
    T <: AbstractDirectProductGroup,
    S₁ <: AbstractAction{<:AbstractGroup{Lie}},
    S₂ <: AbstractAction{<:AbstractGroup{Finite}}
} <: AbstractAction{T}
    group::T
    lie_action::S₁
    finite_action::S₂
end

group(a::DirectProductGroupAction) = a.group
lie_action(a::DirectProductGroupAction) = a.lie_action
finite_action(a::DirectProductGroupAction) = a.finite_action

function Base.show(io::IO, a::DirectProductGroupAction)
    println(io, "DirectProductGroupAction of $(name(group(a)))")
    print(io, " action:")
    for action in actions(a)
        println(io)
        show_action(io, action; offset=2, show_name=true)
    end
end

function act(
    g::GroupElem{T},
    a::MatrixGroupAction{T},
    f::AbstractPolynomialLike
) where T <: AbstractGroup
    sbs = vcat([matrix(g)*v for v in action_vectors(a)]...)
    vars = vcat(action_vectors(a)...)
    return subs(f, vars => sbs)
end

act(
    g::GroupElem{T},
    a::ScalingLieGroupAction{T},
    f::AbstractPolynomialLike
) where T <: AbstractGroup = act(g, MatrixGroupAction(a), f)

function act(
    g::DirectProductGroupElem{T},
    a::AbstractAction{T},
    f::AbstractPolynomialLike
) where T <: AbstractDirectProductGroup
    
end

function are_commutative(
    a₁::AbstractAction{<:AbstractGroup{T₁, F}},
    a₂::AbstractAction{<:AbstractGroup{T₂, F}}
) where {T₁, T₂, F}
    V = space(a₁) + space(a₂)
    v = rand(V)
    g₁, g₂ = rand(group(a₁)), rand(group(a₂))
    v₁ = act(g₂, a₂, act(g₁, a₁, v))
    v₂ = act(g₁, a₁, act(g₂, a₂, v))
    return v₁ ≈ v₂
end

function ⊠(
    a₁::AbstractAction{<:AbstractGroup{T₁, F}},
    a₂::AbstractAction{<:AbstractGroup{T₂, F}}
) where {T₁, T₂, F}
    if are_commutative(a₁, a₂)
        return DirectProductGroupAction(group(a₁) × group(a₂), [a₁, a₂])
    else
        error("Actions are not commutative, cannot form direct product!")
    end
end