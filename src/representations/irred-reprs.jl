export IrreducibleRepresentation


struct IrreducibleRepresentation{
    A<:AbstractGroupAction{Lie}, Wv<:WeightVector
} <: AbstractGroupRepresentation{Lie, HighestWeightModule}
    action::A
    hw_vector::Wv # highest weight vector
end

action(ρ::IrreducibleRepresentation) = ρ.action
hw_vector(ρ::IrreducibleRepresentation) = ρ.hw_vector
space(ρ::IrreducibleRepresentation) = HighestWeightModule(action(ρ), hw_vector(ρ))
highest_weight(ρ::IrreducibleRepresentation) = weight(hw_vector(ρ))
irreducibles(ρ::IrreducibleRepresentation) = [ρ]
isotypic_components(ρ::IrreducibleRepresentation) = [IsotypicComponent(ρ)]

function Base.show(io::IO, ::MIME"text/plain", ρ::IrreducibleRepresentation)
    println(io, "IrreducibleRepresentation of dimension $(dim(ρ))")
    println(io, " Lie group: ", name(group(ρ)))
    print(io, " highest weight: ", highest_weight(ρ))
end

function Base.show(io::IO, ρ::IrreducibleRepresentation)
    print(io, "IrreducibleRepresentation of dimension $(dim(ρ))")
end
