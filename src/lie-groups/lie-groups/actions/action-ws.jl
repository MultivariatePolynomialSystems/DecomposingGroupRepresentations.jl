export hw_spaces, weight_structure

function inv_weight_space(
    a::AbstractGroupAction{Lie, F},
    V::VectorSpace{T, F}
) where {T<:Variable, F}
    @assert variables(a) ⊆ variables(V)
    inv_vars = setdiff(variables(V), variables(a)) # variables invariant under the action
    isempty(inv_vars) && return nothing
    return WeightSpace(zero_weight(group(a)), VariableSpace{T,F}(inv_vars))
end

function weight_structure(
    a::ScalingLieGroupAction{F},
    V::VectorSpace{T, F}
) where {F,S,M,T<:Variable{S,M}}
    @assert variables(a) ⊆ variables(V)
    ws = WeightStructure{F, VectorSpace{T, F}, Weight{weight_type(algebra(a))}}()
    for (i, var) in enumerate(variables(a))
        wv = WeightVector(weight(group(a), i), var)
        push!(ws, wv)
    end
    inv_ws = inv_weight_space(a, V)
    !isnothing(inv_ws) && push!(ws, inv_ws)
    return ws
end

# For the isotypic decomposition V = ⨁ₖ Vₖ returns the weight structure of h.w. spaces of Vⱼ
# If Vₖ = ⨁ⱼ Wⱼ is the irreducible decomposition, then h.w.s.(Vₖ) = ⨁ⱼ h.w.v.(Wⱼ)
hw_spaces(a::ScalingLieGroupAction, V::VectorSpace{<:Variable}) = weight_structure(a, V)


function weight_structure(
    a::MatrixGroupAction{Lie, F},
    V::VectorSpace{<:Variable{T,S}, F};
    as_hw_spaces::Bool=false
) where {F,T,S}
    @assert variables(a) ⊆ variables(V)
    space_type = VectorSpace{Polynomial{T,S,F}, F}
    ws_V = WeightStructure{F, space_type, Weight{weight_type(algebra(a))}}()
    ws_alg = as_hw_spaces ? hw_spaces(algebra(a)) : weight_structure(algebra(a))
    for w_space in ws_alg
        poly_wvs = [sum(av .* vector(wv)) for av in action_vectors(a) for wv in w_space]
        push!(ws_V, WeightSpace(weight(w_space), VectorSpace(F, poly_wvs)))
    end
    inv_ws = inv_weight_space(a, V)
    !isnothing(inv_ws) && push!(ws_V, inv_ws)
    return ws_V
end

hw_spaces(a::MatrixGroupAction{Lie}, V::VectorSpace{<:Variable}) = weight_structure(a, V, as_hw_spaces=true)

function Base.:∩(ws₁::WeightSpace, ws₂::WeightSpace)
    new_weight = vcat(weight(ws₁), weight(ws₂))
    int_space = ∩(space(ws₁), space(ws₂))
    iszero(int_space) && return nothing
    return WeightSpace(new_weight, int_space)
end

function Base.:∩(ws1::WeightStructure, ws2::WeightStructure)
    new_spaces = [w_space₁ ∩ w_space₂ for w_space₁ in ws1 for w_space₂ in ws2]
    return WeightStructure([ws for ws in new_spaces if !isnothing(ws)])
end

hw_spaces(
    a::DirectProductGroupAction{Lie},
    V::VectorSpace{<:Variable}
) = reduce(∩, [hw_spaces(action, V) for action in lie_actions(a)])

function weight_structure(
    a::AbstractGroupAction{Lie, F},
    V::VectorSpace{<:Variable{T,S}, F}
) where {F, T, S}
    hw_struct = hw_spaces(a, V)
    ws = WeightStructure{F, VectorSpace{Polynomial{T,S,F}, F}, Weight{weight_type(algebra(a))}}()
    for hw_space in hw_struct
        for hwv in hw_space
            hwm = HighestWeightModule(a, hwv)
            [push!(ws, wv) for wv in basis(hwm)]
        end
    end
    return ws
end
