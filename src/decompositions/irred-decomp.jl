function ws_nullspace(
    X::AbstractLieAlgebraElem,
    a::AbstractGroupAction{Lie},
    ws::WeightSpace{F};
    tol::Real=1e-5
) where F
    B = basis(space(ws))
    Bₐ = [act(X, a, f) for f in B]
    Cs = zero_combinations(Bₐ)
    isempty(Cs) && return nothing
    return WeightSpace(
        weight(ws),
        VectorSpace(F, [sum(sparsify!(div_by_lowest_magnitude(c, tol), tol) .* B) for c in Cs])
        )
end

function common_nullspace_as_weight_vectors(
    Xs::Vector{<:AbstractLieAlgebraElem},
    a::AbstractGroupAction{Lie},
    ws::WeightSpace
)
    wsₙ = ws
    for X in Xs
        wsₙ = ws_nullspace(X, a, wsₙ)
        isnothing(wsₙ) && return WeightVector[]
    end
    return [WeightVector(weight(wsₙ), v) for v in basis(space(wsₙ))]
end

common_nullspace_as_weight_vectors(
    Xs::Vector{<:AbstractLieAlgebraElem},
    a::AbstractGroupAction{Lie},
    ws::WeightStructure
) = vcat([common_nullspace_as_weight_vectors(Xs, a, w) for w in ws]...)

function irreducibles(
    ρ::GroupRepresentation{A, T}
) where {A<:AbstractGroupAction{Lie}, T<:SymmetricPowersSpace}
    V = space(ρ)
    ws = weight_structure(action(ρ), base_space(V))
    irreds = IrreducibleRepresentation{A}[]
    for p in powers(V)
        sym_ws = sym(ws, p)
        Xs = positive_root_elements(algebra(action(ρ)))
        hw_vectors = common_nullspace_as_weight_vectors(Xs, action(ρ), sym_ws)
        append!(irreds, [IrreducibleRepresentation(action(ρ), hwv) for hwv in hw_vectors])
    end
    return irreds
end

function irreducibles_improved(
    ρ::GroupRepresentation{A, T}
) where {A<:AbstractGroupAction{Lie}, T<:SymmetricPowersSpace}
    V = base_space(space(ρ))
    base_ρ = GroupRepresentation(action(ρ), V)
    irreds = irreducibles(base_ρ)
    
end

function irreducibles(
    ρ::GroupRepresentation{A, T}
) where {A<:AbstractGroupAction{Lie}, T<:AbstractSymmetricPower}
    V = space(ρ)
    ws = weight_structure(action(ρ), base_space(V))
    sym_ws = sym(ws, power(V))
    Xs = positive_root_elements(algebra(action(ρ)))
    hw_vectors = common_nullspace_as_weight_vectors(Xs, action(ρ), sym_ws)
    return [IrreducibleRepresentation(action(ρ), hwv) for hwv in hw_vectors]
end

function irreducibles(
    ρ::GroupRepresentation{A, T}
) where {A<:AbstractGroupAction{Lie}, T<:AbstractDirectSum}
    irreds = IrreducibleRepresentation{A}[]
    for V in summands(space(ρ))
        append!(irreds, irreducibles(GroupRepresentation(action(ρ), V)))
    end
    return irreds
end
