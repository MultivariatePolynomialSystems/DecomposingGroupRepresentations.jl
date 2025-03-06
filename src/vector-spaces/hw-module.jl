export HighestWeightModule


struct HighestWeightModule{T, F, A<:AbstractGroupAction{Lie, F}, Wv<:WeightVector{T}} <: AbstractVectorSpace{F}
    action::A
    hw_vector::Wv
end

action(V::HighestWeightModule) = V.action
group(V::HighestWeightModule) = group(action(V))
hw_vector(V::HighestWeightModule) = V.hw_vector
highest_weight(V::HighestWeightModule) = weight(hw_vector(V))
weight_type(V::HighestWeightModule) = typeof(highest_weight(V))
vector_type(V::HighestWeightModule) = typeof(vector(hw_vector(V)))
variables(V::HighestWeightModule) = free_symbols(vector(hw_vector(V)))
nvariables(V::HighestWeightModule) = length(variables(V))

# Weyl dimension formula
function weyl_dim(λ::Weight, G::AbstractGroup{Lie})
    Δ⁺ = to_vector.(positive_roots(algebra(G)))
    ρ = sum(Δ⁺) / 2
    return Int(prod([dot(λ.weight + ρ, α) / dot(ρ, α) for α in Δ⁺]))
end

dim(V::HighestWeightModule) = weyl_dim(highest_weight(V), group(V))

function Base.show(io::IO, ::MIME"text/plain", V::HighestWeightModule)
    println(io, "HighestWeightModule of dimension ", dim(V))
    println(io, " Lie group: ", name(group(V)))
    print(io, " highest weight: ", highest_weight(V))
end

function Base.show(io::IO, V::HighestWeightModule)
    println(io, "HighestWeightModule of dimension ", dim(V))
end

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

function basis(
    V::HighestWeightModule{T, ComplexF64} where T;
    as_weight_vectors::Bool=false
)
    orb = orbit(hw_vector(V), action(V), Set{weight_type(V)}())
    return as_weight_vectors ? orb : [vector(wv) for wv in orb]
end

function weight_structure(
    V::HighestWeightModule{T, F}
) where {T, F}
    ws = WeightStructure{ExpressionSpace{F}, weight_type(V)}()
    for wv in basis(V; as_weight_vectors=true)
        push!(ws, wv)
    end
    return ws
end
