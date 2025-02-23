export GroupRepresentation,
    irreducibles


struct GroupRepresentation{A<:AbstractGroupAction{Lie}, T<:AbstractVectorSpace} <: AbstractGroupRepresentation{Lie, T}
    action::A
    V::T
end

action(ρ::GroupRepresentation) = ρ.action
group(ρ::GroupRepresentation) = group(action(ρ))
space(ρ::GroupRepresentation) = ρ.V
dim(ρ::GroupRepresentation) = dim(space(ρ))

# called by Shift+Enter
function Base.show(io::IO, ::MIME"text/plain", ρ::GroupRepresentation)
    println(
        io,
        "GroupRepresentation of $(name(group(ρ))) ",
        "on the $(dim(ρ))-dimensional vector space"
    )
    print(io, " Lie group: ", name(group(ρ)))
end

# called by print and inside vectors/matrices
function Base.show(io::IO, ρ::GroupRepresentation)
    print(
        io,
        "GroupRepresentation of $(name(group(ρ))) ",
        "on the $(dim(ρ))-dimensional vector space"
    )
end

function ws_nullspace(X::AbstractLieAlgebraElem, a::AbstractGroupAction{Lie}, ws::WeightSpace{F,T}) where {F,T}
    B = basis(space(ws))
    Bₐ = [act(X, a, f) for f in B]
    Cs = zero_combinations(Bₐ)
    isempty(Cs) && return nothing
    return WeightSpace(weight(ws), T([sum(c .* B) for c in Cs]))
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
) where {A<:AbstractGroupAction{Lie}, T<:AbstractSymmetricPower}
    V = space(ρ)
    ws = weight_structure(action(ρ), base_space(V))
    sym_ws = sym_weight_structure(ws, power(V))
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
