function isotypic_components(
    ρ::GroupRepresentation{A}
) where {A<:AbstractGroupAction{Lie}}
    irreds = irreducibles(ρ)
    W = Weight{weight_type(algebra(action(ρ)))}
    iso_dict = Dict{W, Vector{IrreducibleRepresentation{A}}}()
    for irr in irreds
        irrs = get(iso_dict, highest_weight(irr), nothing)
        if isnothing(irrs)
            iso_dict[highest_weight(irr)] = [irr]
        else
            push!(irrs, irr)
        end
    end
    return [IsotypicComponent(action(ρ), hw, irrs) for (hw, irrs) in iso_dict]
end