export MatrixGroupAction,
    action_vectors,
    ScalingLieGroupAction,
    action_vector,
    DirectProductGroupAction,
    ×,
    are_commutative,
    act

"""
    MatrixGroupAction{T<:GroupType, F} <: AbstractGroupAction{T, F}

Represents a group action of a matrix group on a set of variables.

# Constructors
```julia
MatrixGroupAction(G::S, vectors::AbstractVector{<:AbstractVector{V}}) where {S<:AbstractGroup, V<:Variable}
```

# Examples
```jldoctest
julia> @polyvar x[1:3] y[1:3];

julia> SO3 = LieGroup("SO", 3);

julia> MatrixGroupAction(SO3, [x, y])
MatrixGroupAction of SO(3)
 2 vectors under action: [x₁, x₂, x₃], [y₁, y₂, y₃]
```
"""
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
naction_vectors(a::MatrixGroupAction) = length(action_vectors(a))
DynamicPolynomials.variables(a::MatrixGroupAction) = vcat(action_vectors(a)...)
space(a::MatrixGroupAction{T, F}) where {T,F} = VariableSpace{F}(variables(a))

function Base.show(io::IO, a::MatrixGroupAction)
    println(io, "MatrixGroupAction of $(name(group(a)))")
    print(io, " $(naction_vectors(a)) vectors under action: ")
    print(io, "[", join([join(map(repr, v), ", ") for v in action_vectors(a)], "], ["), "]")
end


# TODO: remove and leave MatrixGroupAction only?
"""
    ScalingLieGroupAction <: AbstractGroupAction

Represents an action of a scaling Lie group on a set of variables.

# Constructors
```julia
ScalingLieGroupAction(v::Vector{<:Variable})
ScalingLieGroupAction(V::AbstractMatrix{<:Variable})
ScalingLieGroupAction(G::ScalingLieGroup, v::Vector{<:Variable})
```

# Examples
```jldoctest
julia> @polyvar x[1:2, 1:3];

julia> ScalingLieGroupAction(x)
ScalingLieGroupAction of (ℂˣ)³
 vector under action: [x₁₋₁, x₂₋₁, x₁₋₂, x₂₋₂, x₁₋₃, x₂₋₃]
 action:
  x₁₋₁ ↦ λ₁x₁₋₁, x₂₋₁ ↦ λ₁x₂₋₁
  x₁₋₂ ↦ λ₂x₁₋₂, x₂₋₂ ↦ λ₂x₂₋₂
  x₁₋₃ ↦ λ₃x₁₋₃, x₂₋₃ ↦ λ₃x₂₋₃

julia> ScalingLieGroupAction(x[:])
ScalingLieGroupAction of ℂˣ
 vector under action: [x₁₋₁, x₂₋₁, x₁₋₂, x₂₋₂, x₁₋₃, x₂₋₃]
 action:
  x₁₋₁ ↦ λx₁₋₁, x₂₋₁ ↦ λx₂₋₁, x₁₋₂ ↦ λx₁₋₂, x₂₋₂ ↦ λx₂₋₂, x₁₋₃ ↦ λx₁₋₃, x₂₋₃ ↦ λx₂₋₃
```
"""
struct ScalingLieGroupAction{F, T <: ScalingLieGroup{F}, V <: Variable} <: AbstractGroupAction{Lie, F}
    group::T
    vars::Vector{V}
end

group(a::ScalingLieGroupAction) = a.group
action_vector(a::ScalingLieGroupAction) = a.vars
DynamicPolynomials.variables(a::ScalingLieGroupAction) = action_vector(a)
space(a::ScalingLieGroupAction{F}) where F = VariableSpace{F}(action_vector(a))

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
    str_action = []
    vars = action_vector(a)
    for j in axes(U, 1)
        nzind, nzval = findnz(U[j, :])
        str_exprs = [repr(λ[j]) * (val != 1 ? superscript(val) : "") * repr(vars[ind]) for (val, ind) in zip(nzval, nzind)]
        push!(str_action, collect(zip(vars[nzind], str_exprs)))
    end
    for free_action in str_action
        print(io, "\n", " "^offset)
        for (j, (var, expr)) in enumerate(free_action)
            print(io, repr(var), " ↦ ", expr)
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

"""
    DirectProductGroupAction <: AbstractGroupAction

Represents an action of a direct product group on a vector space.

# Examples
```jldoctest
julia> @polyvar x[1:3, 1:2];

julia> SO3 = LieGroup("SO", 3);

julia> a₁ = MatrixGroupAction(SO3, eachcol(x))
MatrixGroupAction of SO(3)
 2 vectors under action: [x₁₋₁, x₂₋₁, x₃₋₁], [x₁₋₂, x₂₋₂, x₃₋₂]

julia> a₂ = ScalingLieGroupAction(x)
ScalingLieGroupAction of (ℂˣ)²
 vector under action: [x₁₋₁, x₂₋₁, x₃₋₁, x₁₋₂, x₂₋₂, x₃₋₂]
 action:
  x₁₋₁ ↦ λ₁x₁₋₁, x₂₋₁ ↦ λ₁x₂₋₁, x₃₋₁ ↦ λ₁x₃₋₁
  x₁₋₂ ↦ λ₂x₁₋₂, x₂₋₂ ↦ λ₂x₂₋₂, x₃₋₂ ↦ λ₂x₃₋₂

julia> a₁ × a₂
DirectProductGroupAction of SO(3) × (ℂˣ)²
 lie actions:
```
"""
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