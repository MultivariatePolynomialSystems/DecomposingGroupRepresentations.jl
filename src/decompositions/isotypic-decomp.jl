export IsotypicDecomposition,
    invariant_component,
    representation

"""
    IsotypicDecomposition

Represents the isotypic decomposition of a group representation.

# Constructors
```julia
IsotypicDecomposition(ρ::GroupRepresentation)
```

# Examples
```jldoctest
julia> @polyvar x[1:3] y[1:3];

julia> SO3 = LieGroup("SO", 3);

julia> a₁ = MatrixGroupAction(SO3, [x, y]);

julia> a₂ = ScalingLieGroupAction(vcat(x, y));

julia> a = a₁ × a₂;

julia> V = FixedDegreePolynomials(vcat(x, y), 2);

julia> ρ = GroupRepresentation(a, V)
GroupRepresentation of SO(3) × ℂˣ on 21-dimensional vector space
 Lie group: SO(3) × ℂˣ

julia> isotypics(ρ)
IsotypicDecomposition of SO(3) × ℂˣ-action on 21-dimensional vector space
 number of isotypic components: 3
 multiplicities of irreducibles: 3, 3, 1
 dimensions of irreducibles: 1, 5, 3
 dimensions of isotypic components: 3, 15, 3
```
"""
struct IsotypicDecomposition{A<:AbstractGroupAction{Lie}, T<:GroupRepresentation{A}, W<:Weight, Ic<:IsotypicComponent}
    ρ::T
    iso_dict::Dict{W, Ic}
end

function IsotypicDecomposition(
    irreds::Vector{T},
    ρ::GroupRepresentation{A}
) where {A<:AbstractGroupAction{Lie}, T<:IrreducibleRepresentation{A}}
    a = action(first(irreds))
    W = weight_type(algebra(a))
    iso_dict = Dict{W, IsotypicComponent}()
    for irr in irreds
        hw = highest_weight(irr)
        if haskey(iso_dict, hw)
            push!(iso_dict[hw], irr)
        else
            iso_dict[hw] = IsotypicComponent(a, hw, [irr])
        end
    end
    return IsotypicDecomposition(ρ, iso_dict)
end

function IsotypicDecomposition(
    irr_decomp::IrreducibleDecomposition{A, T, Ir}
) where {A, T, Ir}
    a = action(irr_decomp)
    W = weight_type(algebra(a))
    iso_dict = Dict{W, IsotypicComponent}()
    for irr in irreducibles(irr_decomp)
        hw = highest_weight(irr)
        if haskey(iso_dict, hw)
            push!(iso_dict[hw], irr)
        else
            iso_dict[hw] = IsotypicComponent(a, hw, [irr])
        end
    end
    return IsotypicDecomposition(representation(irr_decomp), iso_dict)
end

function IsotypicDecomposition(
    ρ::GroupRepresentation{A}
) where {A<:AbstractGroupAction{Lie}}
    irr_decomp = irreducibles(ρ)
    return IsotypicDecomposition(irr_decomp)
end

representation(ID::IsotypicDecomposition) = ID.ρ
action(ID::IsotypicDecomposition) = action(ID.ρ)
group(ID::IsotypicDecomposition) = group(ID.ρ)
isotypics(ID::IsotypicDecomposition) = collect(values(ID.iso_dict))
nisotypic(ID::IsotypicDecomposition) = length(ID.iso_dict)
dim(ID::IsotypicDecomposition) = dim(ID.ρ)
Base.getindex(ID::IsotypicDecomposition, w::Weight) = ID.iso_dict[w]
Base.length(ID::IsotypicDecomposition) = length(ID.iso_dict)
Base.iterate(ID::IsotypicDecomposition) = iterate(isotypics(ID))
Base.iterate(ID::IsotypicDecomposition, state) = iterate(isotypics(ID), state)
Base.eltype(::IsotypicDecomposition{A,T,W,Ic}) where {A,T,W,Ic} = Ic

function Base.show(io::IO, iso::IsotypicDecomposition)
    println(io, "IsotypicDecomposition of $(name(group(iso)))-action on $(dim(iso))-dimensional vector space")
    println(io, " number of isotypic components: ", nisotypic(iso))
    if nisotypic(iso) > 50
        println(io, " maximum multiplicity of irreducibles: ", maximum([mul(ic) for ic in isotypics(iso)]))
        print(io, " maximum dimension of isotypic components: ", maximum([dim(ic) for ic in isotypics(iso)]))
    else
        println(io, " multiplicities of irreducibles: ", join([mul(ic) for ic in isotypics(iso)], ", "))
        println(io, " dimensions of irreducibles: ", join([dim_irr(ic) for ic in isotypics(iso)], ", "))
        print(io, " dimensions of isotypic components: ", join([dim(ic) for ic in isotypics(iso)], ", "))
    end
end

invariants(
    ID::IsotypicDecomposition
) = ID[zero(highest_weight(first(isotypics(ID))))]

isotypics(ρ::GroupRepresentation) = IsotypicDecomposition(ρ)