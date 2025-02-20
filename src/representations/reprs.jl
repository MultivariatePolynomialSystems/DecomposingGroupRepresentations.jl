export GroupRepresentation


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