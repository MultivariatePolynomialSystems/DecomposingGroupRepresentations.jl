function isotypic_components(
    irreds::Vector{T}
) where {T<:IrreducibleRepresentation}
    W = weight_type(algebra(action(first(irreds))))
    iso_dict = Dict{W, Vector{T}}()
    for irr in irreds
        irrs = get(iso_dict, highest_weight(irr), nothing)
        if isnothing(irrs)
            iso_dict[highest_weight(irr)] = [irr]
        else
            push!(irrs, irr)
        end
    end
    return [IsotypicComponent(action(first(irreds)), hw, irrs) for (hw, irrs) in iso_dict]
end

function isotypic_components(
    ρ::GroupRepresentation{A}
) where {A<:AbstractGroupAction{Lie}}
    irreds = irreducibles(ρ)
    return isotypic_components(irreds)
end