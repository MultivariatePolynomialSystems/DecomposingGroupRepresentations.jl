export HighestWeightModule


struct HighestWeightModule{F, A<:AbstractGroupAction{Lie, F}, Wv<:WeightVector} <: AbstractVectorSpace{F}
    action::A
    hw_vector::Wv
end

action(V::HighestWeightModule) = V.action
group(V::HighestWeightModule) = group(action(V))
hw_vector(V::HighestWeightModule) = V.hw_vector
highest_weight(V::HighestWeightModule) = weight(hw_vector(V))
weight_type(V::HighestWeightModule) = typeof(highest_weight(V))

# Weyl dimension formula
function weyl_dim(λ::Weight, G::AbstractGroup{Lie})
    Δ⁺ = Vector{Int}.(positive_roots(algebra(G)))
    ρ = sum(Δ⁺) / 2
    return Int(prod([dot(λ + ρ, α) / dot(ρ, α) for α in Δ⁺]))
end

dim(V::HighestWeightModule) = weyl_dim(highest_weight(V), group(V))

function orbit(
    wv::T,
    action::AbstractGroupAction{Lie},
    processed_weights::Set{<:Weight};
    tol::Real=1e-5
) where T<:WeightVector
    isapprox(vector(wv), zero(vector(wv)); atol=tol) && return T[]
    push!(processed_weights, weight(wv))
    orbit_vectors = [wv]
    for nre in negative_root_elements(algebra(action))
        if weight(wv) + root(nre) ∉ processed_weights
            new_wv = act(nre, action, wv)
            append!(orbit_vectors, orbit(new_wv, action, processed_weights; tol=tol))
        end
    end
    return orbit_vectors
end

basis(V::HighestWeightModule{ComplexF64}) = orbit(hw_vector(V), action(V), Set{weight_type(V)}())