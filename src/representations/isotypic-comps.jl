export IsotypicComponent


struct IsotypicComponent{
    A<:AbstractGroupAction{Lie}, W<:Weight, Ir<:IrreducibleRepresentation{A}
} <: AbstractGroupRepresentation{Lie, DirectSumSpace}
    action::A
    highest_weight::W
    irreds::Vector{Ir}
end

action(ρ::IsotypicComponent) = ρ.action
group(ρ::IsotypicComponent) = group(action(ρ))
highest_weight(ρ::IsotypicComponent) = ρ.highest_weight
irreducibles(ρ::IsotypicComponent) = ρ.irreds
nirreducible(ρ::IsotypicComponent) = length(ρ.irreds)
space(ρ::IsotypicComponent) = DirectSumSpace([space(ρᵢ) for ρᵢ in irreducibles(ρ)])
dim(ρ::IsotypicComponent) = sum(dim(ρᵢ) for ρᵢ in irreducibles(ρ))
hw_vectors(ρ::IsotypicComponent) = [hw_vector(ρᵢ) for ρᵢ in irreducibles(ρ)]

function Base.show(io::IO, ::MIME"text/plain", ρ::IsotypicComponent)
    println(io, "IsotypicComponent of dimension ", dim(ρ))
    println(io, " Lie group: ", name(group(ρ)))
    println(io, " highest weight: ", highest_weight(ρ))
    println(io, " multiplicity of irreducible subrepresentation: ", nirreducible(ρ))
    print(io, " dimension of irreducible subrepresentation: ", weyl_dim(highest_weight(ρ)))
end

function Base.show(io::IO, ρ::IsotypicComponent)
    print(io, "IsotypicComponent of dimension ", dim(ρ))
end