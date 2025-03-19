export irreducibles, irreducibles_old

function ws_nullspace(
    X::AbstractLieAlgebraElem,
    a::AbstractGroupAction{Lie},
    ws::WeightSpace;
    tol::Real=1e-5
)
    B = basis(space(ws))
    Bₐ = [act(X, a, f) for f in B]
    if all(f -> isapprox(f, zero(f), atol=tol), Bₐ)
        return ws
    end
    dim(ws) == 1 && return nothing
    Cs = zero_combinations(Bₐ)
    isempty(Cs) && return nothing
    return WeightSpace(
        weight(ws),
        VectorSpace(field_type(ws), [sum(c .* B) for c in Cs]; in_rref=false)
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
) where {A<:AbstractGroupAction{Lie}, T<:DirectSum}
    irreds = IrreducibleRepresentation{A}[]
    for V in summands(space(ρ))
        append!(irreds, irreducibles(GroupRepresentation(action(ρ), V)))
    end
    return irreds
end

function irreducibles(
    ρ::GroupRepresentation{A, S}
) where {A<:AbstractGroupAction{Lie}, S<:VectorSpace{<:Variable}}
    ws = weight_structure(action(ρ), space(ρ))
    Xs = positive_root_elements(algebra(action(ρ)))
    hw_vectors = common_nullspace_as_weight_vectors(Xs, action(ρ), ws)
    return [IrreducibleRepresentation(action(ρ), hwv) for hwv in hw_vectors]
end

function irreducibles(
    ρ::GroupRepresentation{A, S}
) where {T, F, A<:AbstractGroupAction{Lie}, S<:SymmetricPower{T, F, <:HighestWeightModule}}
    V = space(ρ)
    Vb = base_space(V)
    power(V) == 1 && return [IrreducibleRepresentation(action(ρ), Vb)]
    ws = weight_structure(Vb)
    sym_ws = sym_weight_struct(highest_weight(Vb), algebra(action(ρ)), power(V), ws)
    println("new nweights: ", nweights(sym_ws))
    println("new sum of dims: ", sum([dim(space(sym_ws[w])) for w in weights(sym_ws)]))
    Xs = positive_root_elements(algebra(action(ρ)))
    hw_vectors = common_nullspace_as_weight_vectors(Xs, action(ρ), sym_ws)
    return [IrreducibleRepresentation(action(ρ), hwv) for hwv in hw_vectors]
end

function irreducibles_old(
    ρ::GroupRepresentation{A, S}
) where {T, F, A<:AbstractGroupAction{Lie}, S<:SymmetricPower{T, F, <:HighestWeightModule}}
    V = space(ρ)
    power(V) == 1 && return [IrreducibleRepresentation(action(ρ), base_space(V))]
    ws = weight_structure(base_space(V))
    sym_ws = sym(ws, power(V))
    println("old nweights: ", nweights(sym_ws))
    println("old sum of dims: ", sum([dim(space(sym_ws[w])) for w in weights(sym_ws)]))
    Xs = positive_root_elements(algebra(action(ρ)))
    hw_vectors = common_nullspace_as_weight_vectors(Xs, action(ρ), sym_ws)
    return [IrreducibleRepresentation(action(ρ), hwv) for hwv in hw_vectors]
end

function irreducibles_old(
    ρ::GroupRepresentation{A, T}
) where {A<:AbstractGroupAction{Lie}, T<:SymmetricPower}
    V = space(ρ)
    ρ_base = GroupRepresentation(action(ρ), base_space(V))
    irreds_base = irreducibles(ρ_base)
    if length(irreds_base) == 1
        V = SymmetricPower(space(first(irreds_base)), power(V))
        ρᵢ = GroupRepresentation(action(ρ), V)
        return irreducibles_old(ρᵢ)
    end
    mexps = multiexponents(degree=power(V), nvars=length(irreds_base))
    irreds = IrreducibleRepresentation{A}[]
    for mexp in mexps
        Vᵢ = TensorProduct([SymmetricPower(space(irreds_base[idx]), val) for (val, idx) in zip(mexp.nzval, mexp.nzind)])
        ρᵢ = GroupRepresentation(action(ρ), Vᵢ)
        append!(irreds, irreducibles(ρᵢ))
    end
    return irreds
end

function irreducibles(
    ρ::GroupRepresentation{A, T}
) where {A<:AbstractGroupAction{Lie}, T<:SymmetricPower}
    V = space(ρ)
    ρ_base = GroupRepresentation(action(ρ), base_space(V))
    irreds_base = irreducibles(ρ_base)
    if length(irreds_base) == 1
        V = SymmetricPower(space(first(irreds_base)), power(V))
        ρᵢ = GroupRepresentation(action(ρ), V)
        return irreducibles(ρᵢ)
    end
    mexps = multiexponents(degree=power(V), nvars=length(irreds_base))
    irreds = IrreducibleRepresentation{A}[]
    for mexp in mexps
        Vᵢ = TensorProduct([SymmetricPower(space(irreds_base[idx]), val) for (val, idx) in zip(mexp.nzval, mexp.nzind)])
        ρᵢ = GroupRepresentation(action(ρ), Vᵢ)
        append!(irreds, irreducibles(ρᵢ))
    end
    return irreds
end

function irreducibles(
    ρ::GroupRepresentation{A, S}
) where {T, F, A<:AbstractGroupAction{Lie}, S<:TensorProduct{T, F, <:HighestWeightModule}}
    V = space(ρ)
    if nfactors(V) == 2
        ws₁, ws₂ = weight_structure(factor(V, 1)), weight_structure(factor(V, 2)) # TODO: save once computed
        println("nweights 1: ", nweights(ws₁))
        println("nweights 2: ", nweights(ws₂))
        tensor_ws = tensor(ws₁, ws₂)
        println("nweights 1 ⊗ 2: ", nweights(tensor_ws))
        Xs = positive_root_elements(algebra(action(ρ)))
        hw_vectors = common_nullspace_as_weight_vectors(Xs, action(ρ), tensor_ws)
        return [IrreducibleRepresentation(action(ρ), hwv) for hwv in hw_vectors]
    end
    W = TensorProduct(factors(V, [1,2]))
    ρW = GroupRepresentation(action(ρ), W)
    irrs = irreducibles(ρW)
    irreds = IrreducibleRepresentation{A}[]
    for irr in irrs
        Vᵣ = TensorProduct(space(irr), factors(V, 3:nfactors(V)))
        ρᵣ = GroupRepresentation(action(ρ), Vᵣ)
        append!(irreds, irreducibles(ρᵣ))
    end
    return irreds
end

function irreducibles(
    ρ::GroupRepresentation{A, S}
) where {A<:AbstractGroupAction{Lie}, S<:TensorProduct}
    if nfactors(space(ρ)) == 2
        irrs₁ = irreducibles(GroupRepresentation(action(ρ), factor(space(ρ), 1)))
        irrs₂ = irreducibles(GroupRepresentation(action(ρ), factor(space(ρ), 2)))
        irreds = IrreducibleRepresentation{A}[]
        for irr₁ in irrs₁, irr₂ in irrs₂
            V = TensorProduct([space(irr₁), space(irr₂)])
            ρV = GroupRepresentation(action(ρ), V)
            append!(irreds, irreducibles(ρV))
        end
        return irreds
    end
    V = TensorProduct(factors(space(ρ), 2:nfactors(space(ρ))))
    ρV = GroupRepresentation(action(ρ), V)
    irrs = irreducibles(ρV)
    irreds = IrreducibleRepresentation{A}[]
    for irr in irrs
        ρ₁ = GroupRepresentation(action(ρ), factor(space(ρ), 1))
        irrs₁ = irreducibles(ρ₁)
        for irr₁ in irrs₁
            W = TensorProduct([space(irr₁), space(irr)])
            ρW = GroupRepresentation(action(ρ), W)
            append!(irreds, irreducibles(ρW))
        end
    end
    return irreds
end
