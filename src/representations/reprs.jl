export GroupRepresentation,
    irreducibles,
    common_nullspace_as_weight_vectors, ws_nullspace,
    isotypic_components


struct GroupRepresentation{A<:AbstractGroupAction{Lie}, T<:AbstractSpace} <: AbstractGroupRepresentation{Lie, T}
    action::A
    V::T
end

action(ρ::GroupRepresentation) = ρ.action
space(ρ::GroupRepresentation) = ρ.V

# called by Shift+Enter
function Base.show(io::IO, ::MIME"text/plain", ρ::GroupRepresentation)
    println(
        io,
        "GroupRepresentation of $(name(group(ρ))) ",
        "on $(dim(ρ))-dimensional vector space"
    )
    print(io, " Lie group: ", name(group(ρ)))
end

# called by print and inside vectors/matrices
function Base.show(io::IO, ρ::GroupRepresentation)
    print(
        io,
        "GroupRepresentation of $(name(group(ρ))) ",
        "on $(dim(ρ))-dimensional vector space"
    )
end

summands(
    ρ::GroupRepresentation{<:AbstractGroupAction, <:AbstractDirectSum}
) = [GroupRepresentation(action(ρ), V) for V in summands(space(ρ))]

GroupRepresentation(a::AbstractGroupAction{Lie}, V::HighestWeightModule) = IrreducibleRepresentation(a, V)