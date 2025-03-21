export IrreducibleRepresentation,
    highest_weight,
    hw_vector


struct IrreducibleRepresentation{
    A<:AbstractGroupAction{Lie}, T<:HighestWeightModule
} <: AbstractGroupRepresentation{Lie, T}
    action::A
    hw_module::T
end

IrreducibleRepresentation(
    action::AbstractGroupAction{Lie},
    hw_vector::WeightVector
) = IrreducibleRepresentation(action, HighestWeightModule(action, hw_vector))

action(ρ::IrreducibleRepresentation) = ρ.action
hw_vector(ρ::IrreducibleRepresentation) = hw_vector(ρ.hw_module)
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
