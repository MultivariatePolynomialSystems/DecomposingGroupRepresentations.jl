export AbstractLieAlgebraAction,
    LieAlgebraAction,
    ScalingLieAction,
    SumLieAlgebraAction,
    var_groups,
    hw_spaces


act(
    X::LieAlgebraElem{F, T},
    a::MatrixGroupAction{Lie, F},
    f::AbstractPolynomialLike
) where {F, T<:AbstractLieAlgebra{F}} = expand(sum([dot(differentiate(f, vars), -matrix(X)*vars) for vars in action_vectors(a)]))

act(
    X::LieAlgebraElem{F, ScalingLieAlgebra{F}},
    a::ScalingLieGroupAction{F},
    f::AbstractPolynomialLike
) where {F} = act(X, MatrixGroupAction(a), f)

function weight_structure(a::LieAlgebraAction, vars::Vector{Variable}; as_hw_spaces::Bool=false)
    @assert variables(a) ⊆ vars
    vars_dict = Dict(zip(vars, 1:length(vars)))
    ws = WeightStructure()
    alg_ws = as_hw_spaces ? hw_spaces(algebra(a)) : weight_structure(algebra(a))
    for w_space in alg_ws
        new_w_space = WeightSpace(weight(w_space), zeros(ComplexF64, length(vars), dim(w_space)*nvar_groups(a)))
        for (i, vars_group) in enumerate(var_groups(a))
            for (j, var) in enumerate(vars_group)
                new_w_space[vars_dict[var], (i-1)*dim(w_space)+1:i*dim(w_space)] = w_space[j, :]
            end
        end
        push!(ws, new_w_space)
    end
    inv_ws = inv_weight_space(a, vars)
    !isnothing(inv_ws) && push!(ws, inv_ws)
    return ws
end

hw_spaces(a::LieAlgebraAction, vars::Vector{Variable}) = weight_structure(a, vars, as_hw_spaces=true)


struct SumLieAlgebraAction <: AbstractLieAlgebraAction
    alg::SumLieAlgebra
    actions::Vector{AbstractLieAlgebraAction}
end

algebra(a::SumLieAlgebraAction) = a.alg
actions(a::SumLieAlgebraAction) = a.actions
nsummands(a::SumLieAlgebraAction) = length(actions(a))
Base.getindex(a::SumLieAlgebraAction, i::Int) = actions(a)[i]

function Base.show(io::IO, a::SumLieAlgebraAction)
    println(io, "SumLieAlgebraAction of $(name(algebra(a)))")
    print(io, " action:")
    for action in actions(a)
        println(io)
        show_action(io, action; offset=2, show_name=true)
    end
end

# TODO
function are_commutative(a₁::AbstractLieAlgebraAction, a₂::AbstractLieAlgebraAction)
    return true
end

function ⊕(a₁::AbstractLieAlgebraAction, a₂::AbstractLieAlgebraAction)
    @assert are_commutative(a₁, a₂)
    alg = algebra(a₁) ⊕ algebra(a₂)
    return SumLieAlgebraAction(alg, [a₁, a₂])
end

function ⊕(a₁::SumLieAlgebraAction, a₂::AbstractLieAlgebraAction)
    @assert are_commutative(a₁, a₂)
    alg = algebra(a₁) ⊕ algebra(a₂)
    return SumLieAlgebraAction(alg, [actions(a₁)..., a₂])
end

function act(elem::SumLieAlgebraElem, f::Union{Expression, Monomial}, action::SumLieAlgebraAction)
    return expand(sum([act(elem[i], f, action[i]) for i in 1:nsummands(action)]))
end

hw_spaces(
    a::SumLieAlgebraAction,
    vars::Vector{Variable}
) = ∩([hw_spaces(action, vars) for action in actions(a)])

function weight_structure(a::AbstractLieAlgebraAction, vars::Vector{Variable})
    hws = hw_spaces(a, vars)
    mons = MonomialBasis{Int8,Int16}(vars, degree=1, upto=false)
    ws = WeightStructure()
    for hw_space in hws
        for hwv in hw_space
            wexpr = WeightExpression(hwv, mons)
            orb = orbit(wexpr, a)
            [push!(ws, WeightVector(we)) for we in orb]
        end
    end
    return ws
end

function as_matrix(elem::AbstractLieAlgebraElem, B::MonomialBasis, action::AbstractLieAlgebraAction)
    @assert algebra(elem) == algebra(action)
    M = zeros(ComplexF64, length(B), length(B))
    for (i, mon) in enumerate(B)
        gMon = act(elem, mon, action)
        M[:, i] = coefficients(gMon, B)
    end
    return M
end