export IsotypicDecomposition,
    invariant_component

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
    irr_decomp = IrreducibleDecomposition(ρ)
    return IsotypicDecomposition(irr_decomp)
end

representation(ID::IsotypicDecomposition) = ID.ρ
action(ID::IsotypicDecomposition) = action(ID.ρ)
group(ID::IsotypicDecomposition) = group(ID.ρ)
isotypic_components(ID::IsotypicDecomposition) = values(ID.iso_dict)
nisotypic(ID::IsotypicDecomposition) = length(ID.iso_dict)
dim(ID::IsotypicDecomposition) = dim(ID.ρ)
Base.getindex(ID::IsotypicDecomposition, w::Weight) = ID.iso_dict[w]
Base.length(ID::IsotypicDecomposition) = length(ID.iso_dict)
Base.iterate(ID::IsotypicDecomposition) = iterate(isotypic_components(ID))
Base.iterate(ID::IsotypicDecomposition, state) = iterate(isotypic_components(ID), state)

function Base.show(io::IO, iso::IsotypicDecomposition)
    println(io, "IsotypicDecomposition of $(name(group(iso)))-action on $(dim(iso))-dimensional vector space")
    println(io, " number of isotypic components: ", nisotypic(iso))
    println(io, " multiplicities of irreducibles: ", join([mul(ic) for ic in isotypic_components(iso)], ", "))
    print(io, " dimensions of isotypic components: ", join([dim(ic) for ic in isotypic_components(iso)], ", "))
end

invariant_component(
    ID::IsotypicDecomposition
) = ID[zero(highest_weight(first(isotypic_components(ID))))]
