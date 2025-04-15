export IrreducibleRepresentation,
    highest_weight,
    hw_vector


"""
    IrreducibleRepresentation <: AbstractGroupRepresentation

Represents an irreducible group representation.

# Examples
```jldoctest
julia> @polyvar x y z;

julia> vars = [x, y, z];

julia> SO3 = LieGroup("SO", 3);

julia> a = MatrixGroupAction(SO3, [vars]);

julia> V = FixedDegreePolynomials(vars, 2);

julia> ρ = GroupRepresentation(a, V);

julia> irr_decomp = irreducibles(ρ)
IrreducibleDecomposition of SO(3)-action on 6-dimensional vector space
 number of irreducibles: 2
 dimensions of irreducibles: 1, 5

julia> irr_decomp[2]
IrreducibleRepresentation of dimension 5
 Lie group: SO(3)
 highest weight: [2]

julia> space(irr_decomp[2])
HighestWeightModule of dimension 5
 Lie group: SO(3)
 highest weight: [2]

julia> basis(space(irr_decomp[2]))
5-element Vector{Polynomial{Commutative{CreationOrder}, Graded{LexOrder}, ComplexF64}}:
 -y² + (0.0 + 2.0im)xy + x²
 (0.0 + 1.0im)yz + xz
 (2.0 + 0.0im)z² + -y² + -x²
 (0.0 + 1.0im)yz + -xz
 -y² + (-0.0 - 2.0im)xy + x²
```
"""
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
irreducibles(ρ::IrreducibleRepresentation; decomp_count=nothing) = [ρ]
isotypics(ρ::IrreducibleRepresentation) = [IsotypicComponent(ρ)]

function Base.show(io::IO, ::MIME"text/plain", ρ::IrreducibleRepresentation)
    println(io, "IrreducibleRepresentation of dimension $(dim(ρ))")
    println(io, " Lie group: ", name(group(ρ)))
    print(io, " highest weight: ", highest_weight(ρ))
end

function Base.show(io::IO, ρ::IrreducibleRepresentation)
    print(io, "IrreducibleRepresentation of dimension $(dim(ρ))")
end
