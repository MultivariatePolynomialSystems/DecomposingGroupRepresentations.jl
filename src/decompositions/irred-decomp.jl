export irreducibles, irreducibles_old, IrreducibleDecomposition

function ws_nullspace(
    X::AbstractLieAlgebraElem,
    a::AbstractGroupAction{Lie},
    ws::WeightSpace;
    tol::Real=1e-5,
    logging::Bool=true
)
    B = basis(space(ws))
    Bₐ = [act(X, a, f) for f in B]

    logging && print(crayon"#f4d03f", "Applied 'act' method $(dim(space(ws))) times\n")

    if all(f -> isapprox(f, zero(f), atol=tol), Bₐ)
        return ws
    end
    dim(ws) == 1 && return nothing
    Cs = zero_combinations(Bₐ)
    isempty(Cs) && return nothing
    return WeightSpace(
        weight(ws),
        PolySpace{field_type(ws)}([div_by_smallest_coeff(sum(c .* B)) for c in Cs])
        )
end

function common_nullspace_as_weight_vectors(
    Xs::Vector{<:AbstractLieAlgebraElem},
    a::AbstractGroupAction{Lie},
    ws::WeightSpace,
    nhwv::Int;
    logging::Bool=true
)
    dim(ws) == 1 && return [WeightVector(weight(ws), v) for v in basis(space(ws))]

    logging && printstyled("Looking for the highest weight vectors in $(str_irr("W", weight(ws))) of dimension $(dim(ws))...\n", color=:green)
    wsₙ = ws
    for X in Xs
        wsₙ = ws_nullspace(X, a, wsₙ)
        isnothing(wsₙ) && return WeightVector[]
    end

    @assert dim(wsₙ) == nhwv
    logging && printstyled("Found: ", color=:green)
    logging && printstyled("$(join(map(repr, basis(space(wsₙ))), ", "))\n", color=:green)

    return [WeightVector(weight(wsₙ), div_by_smallest_coeff(v)) for v in basis(space(wsₙ))] # TODO: remove div_by_smallest_coeff
end

function common_nullspace_as_weight_vectors(
    Xs::Vector{<:AbstractLieAlgebraElem},
    a::AbstractGroupAction{Lie},
    ws::WeightStructure,
    cg_decomp::Dict{<:Weight, Int}
)
    return vcat([begin
        common_nullspace_as_weight_vectors(Xs, a, ws[w], cg_decomp[w])
    end for w in weights(ws)]...)
end

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
    ρ::GroupRepresentation{A, S};
    logging::Bool=true
) where {A<:AbstractGroupAction{Lie}, S<:VariableSpace}
    logging && printstyled("Decomposing the space of variables $(to_string(space(ρ))) into irreducibles...\n", color=:red)

    hws = hw_spaces(action(ρ), space(ρ))
    hw_vectors = [hwv for ws in hws for hwv in ws]

    logging && printstyled("⟨$(join(map(repr, basis(space(ρ))), ", "))⟩ = ", color=:red)
    logging && printstyled(join([str_irr("V", weight(hwv)) for hwv in hw_vectors], " ⊕ "), "\n", color=:red)
    logging && printstyled("$(length(hw_vectors)) irreducibles\n", color=:red)
    
    return [IrreducibleRepresentation(action(ρ), hwv) for hwv in hw_vectors]
end

function irreducibles(
    ρ::GroupRepresentation{A, S};
    logging::Bool=true
) where {T, F, A<:AbstractGroupAction{Lie}, S<:SymmetricPower{T, F, <:HighestWeightModule}}
    V = space(ρ)
    logging && printstyled("Decomposing $(to_string(V)) into irreducibles...\n", color=:blue)

    Vb = base_space(V)
    power(V) == 1 && return [IrreducibleRepresentation(action(ρ), Vb)]
    ws = weight_structure(Vb)
    
    logging && printstyled("Applying Clebsch-Gordan to decompose $(to_string(V)) into irreducibles...\n", color=:blue)
    sym_ws, cg_decomp = sym_weight_struct(highest_weight(Vb), algebra(action(ρ)), power(V), ws)
    logging && printstyled(to_string(V), " = ", str_irr_decomp("V", cg_decomp), "\n", color=:blue)

    Xs = positive_root_elements(algebra(action(ρ)))
    hw_vectors = common_nullspace_as_weight_vectors(Xs, action(ρ), sym_ws, cg_decomp)
    return [IrreducibleRepresentation(action(ρ), hwv) for hwv in hw_vectors]
end

function irreducibles(
    ρ::GroupRepresentation{A, T};
    logging::Bool=true
) where {A<:AbstractGroupAction{Lie}, T<:SymmetricPower}
    V = space(ρ)
    logging && printstyled("Decomposing $(to_string(V)) into irreducibles...\n", color=:magenta)

    ρ_base = GroupRepresentation(action(ρ), base_space(V))
    irreds_base = irreducibles(ρ_base)
    if length(irreds_base) == 1
        logging && printstyled("$(to_string(base_space(V))) is already irreducible, decomposing its symmetric power directly...\n", color=:magenta)
        V = SymmetricPower(space(first(irreds_base)), power(V))
        ρᵢ = GroupRepresentation(action(ρ), V)
        return irreducibles(ρᵢ)
    end
    mexps = multiexponents(degree=power(V), nvars=length(irreds_base))

    Vs = [TensorProduct([SymmetricPower(space(irreds_base[idx]), val) for (val, idx) in zip(mexp.nzval, mexp.nzind)]) for mexp in mexps]
    if logging
        printstyled("Sym$(superscript(power(V)))(", color=:magenta)
        printstyled(join([to_string(space(irr)) for irr in irreds_base], " ⊕ "), ") = ", color=:magenta)
        printstyled(join(["(" * to_string(V) * ")" for V in Vs], " ⊕ "), color=:magenta)
        println()
    end

    irreds = IrreducibleRepresentation{A}[]
    for Vᵢ in Vs
        ρᵢ = GroupRepresentation(action(ρ), Vᵢ)
        append!(irreds, irreducibles(ρᵢ))
    end
    return irreds
end

function irreducibles(
    ρ::GroupRepresentation{A, S};
    logging::Bool=true
) where {T, F, A<:AbstractGroupAction{Lie}, S<:TensorProduct{T, F, <:HighestWeightModule, <:HighestWeightModule}}
    logging && printstyled("Decomposing $(to_string(V)) into irreducibles...\n", color=:blue)
    V = space(ρ)
    hw₁, hw₂ = highest_weight(first(V)), highest_weight(second(V))
    
    logging && printstyled("Applying Clebsch-Gordan to decompose $(to_string(V)) into irreducibles...\n", color=:blue)
    ws₁, ws₂ = weight_structure(first(V)), weight_structure(second(V)) # TODO: save once computed
    tensor_ws, cg_decomp = tensor_weight_struct(hw₁, hw₂, algebra(action(ρ)), ws₁, ws₂)
    logging && printstyled(to_string(V), " = ", join([str_irr("V", w) for w in weights(tensor_ws)], " ⊕ "), "\n", color=:blue)

    Xs = positive_root_elements(algebra(action(ρ)))
    hw_vectors = common_nullspace_as_weight_vectors(Xs, action(ρ), tensor_ws, cg_decomp)
    return [IrreducibleRepresentation(action(ρ), hwv) for hwv in hw_vectors]
end

function irreducibles(
    ρ::GroupRepresentation{A, S};
    logging::Bool=true
) where {A<:AbstractGroupAction{Lie}, S<:TensorProduct}
    V = space(ρ)
    logging && printstyled("Decomposing $(to_string(V)) into irreducibles...\n", color=:magenta)
    irrs₁ = irreducibles(GroupRepresentation(action(ρ), first(space(ρ))))
    irrs₂ = irreducibles(GroupRepresentation(action(ρ), second(space(ρ))))

    logging && printstyled("Decomposing tensor product knowing irreducibles of the factors...\n", color=:magenta)
    logging && printstyled("($(join([to_string(space(irr)) for irr in irrs₁], " ⊕ "))) ⊗ ", color=:magenta)
    logging && printstyled("($(join([to_string(space(irr)) for irr in irrs₂], " ⊕ "))) = ", color=:magenta)
    logging && printstyled(join(["($(to_string(space(irr₁))) ⊗ $(to_string(space(irr₂))))" for (irr₁, irr₂) in product(irrs₁, irrs₂)], " ⊕ "), "\n", color=:magenta)

    irreds = IrreducibleRepresentation{A}[]
    for irr₁ in irrs₁, irr₂ in irrs₂
        V = TensorProduct([space(irr₁), space(irr₂)])
        ρV = GroupRepresentation(action(ρ), V)
        append!(irreds, irreducibles(ρV))
    end
    return irreds
end

function irreducibles(
    ρ::GroupRepresentation{A, S}
) where {A<:AbstractGroupAction{Lie}, S<:FixedDegreePolynomials}
    V = space(ρ)
    return irreducibles(GroupRepresentation(action(ρ), SymmetricPower(var_space(V), degree(V))))
end

function irreducibles(
    ρ::GroupRepresentation{A, S}
) where {A<:AbstractGroupAction{Lie}, S<:FixedMultidegreePolynomials}
    V = space(ρ)
    return irreducibles(GroupRepresentation(action(ρ), TensorProduct([SymmetricPower(Vᵢ, dᵢ) for (Vᵢ, dᵢ) in zip(var_spaces(V), degrees(V))])))
end


struct IrreducibleDecomposition{A<:AbstractGroupAction{Lie}, T<:GroupRepresentation{A}, Ir<:IrreducibleRepresentation{A}}
    ρ::T
    irreds::Vector{Ir}

    function IrreducibleDecomposition(
        ρ::T,
        irreds::Vector{Ir}
    ) where {A<:AbstractGroupAction{Lie}, T<:GroupRepresentation{A}, Ir<:IrreducibleRepresentation{A}}
        @assert dim(space(ρ)) == sum([dim(space(irr)) for irr in irreds])
        new{A, T, Ir}(ρ, irreds)
    end
end

IrreducibleDecomposition(
    ρ::GroupRepresentation
) = IrreducibleDecomposition(ρ, irreducibles(ρ))

representation(ID::IrreducibleDecomposition) = ID.ρ
action(ID::IrreducibleDecomposition) = action(ID.ρ)
group(ID::IrreducibleDecomposition) = group(ID.ρ)
irreducibles(ID::IrreducibleDecomposition) = ID.irreds
nirreducibles(ID::IrreducibleDecomposition) = length(irreducibles(ID))
dim(ID::IrreducibleDecomposition) = dim(ID.ρ)
Base.getindex(ID::IrreducibleDecomposition, i::Int) = ID.irreds[i]

function Base.show(io::IO, ID::IrreducibleDecomposition)
    println(io, "IrreducibleDecomposition of $(name(group(ID)))-action on $(dim(ID))-dimensional vector space")
    println(io, " number of irreducibles: ", nirreducibles(ID))
    print(io, " dimensions of irreducibles: ", join([dim(irr) for irr in irreducibles(ID)], ", "))
end