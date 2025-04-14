export HighestWeightModule,
    weyl_dim


mutable struct HighestWeightModule{T, F, A<:AbstractGroupAction{Lie, F}, Wv<:WeightVector{T}} <: AbstractSpace{T, F}
    action::A
    hw_vector::Wv
    # basis::Union{Vector{Wv}, Nothing}
end

# HighestWeightModule(
#     action::A,
#     hw_vector::Wv
# ) where {A<:AbstractGroupAction{Lie}, Wv<:WeightVector} = HighestWeightModule(action, hw_vector, nothing)

action(V::HighestWeightModule) = V.action
group(V::HighestWeightModule) = group(action(V))
algebra(V::HighestWeightModule) = algebra(action(V))
hw_vector(V::HighestWeightModule) = V.hw_vector
highest_weight(V::HighestWeightModule) = weight(hw_vector(V))
weight_type(V::HighestWeightModule) = typeof(highest_weight(V))
vector_type(V::HighestWeightModule) = typeof(vector(hw_vector(V)))
DynamicPolynomials.variables(V::HighestWeightModule) = variables(vector(hw_vector(V)))
DynamicPolynomials.nvariables(V::HighestWeightModule) = nvariables(vector(hw_vector(V)))

# Weyl dimension formula
function weyl_dim(λ::Weight, A::AbstractLieAlgebra)
    Δ⁺ = to_vector.(positive_roots(A))
    ρ = sum(Δ⁺) / 2
    return Int(prod([dot(λ.weight + ρ, α) / dot(ρ, α) for α in Δ⁺]))
end

dim(V::HighestWeightModule) = weyl_dim(highest_weight(V), algebra(V))

function Base.show(io::IO, ::MIME"text/plain", V::HighestWeightModule; indent::Int=0)
    println(io, " "^indent, "HighestWeightModule of dimension ", dim(V))
    println(io, " "^indent, " Lie group: ", name(group(V)))
    print(io, " "^indent, " highest weight: ", highest_weight(V))
end

function Base.show(io::IO, V::HighestWeightModule; indent::Int=0)
    print(io, " "^indent, "HighestWeightModule of dimension ", dim(V))
end

function orbit(
    wv::T,
    action::AbstractGroupAction{Lie},
    processed_weights::Set{<:Weight};
    tol::Real=1e-5
) where T<:WeightVector
    isapprox(vector(wv), zero(vector(wv)); atol=tol) && return T[]
    wv = WeightVector(weight(wv), div_by_smallest_coeff(vector(wv)))
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
    # if !isnothing(V.basis)
    #     B = V.basis
    # else
    #     B = orbit(hw_vector(V), action(V), Set{weight_type(V)}())
    #     V.basis = B
    # end
    B = orbit(hw_vector(V), action(V), Set{weight_type(V)}())
    return as_weight_vectors ? B : [vector(wv) for wv in B]
end

function weight_structure(
    V::HighestWeightModule{T, F}
) where {T<:Polynomial, F}
    ws = WeightStructure{PolySpace{F, T}, weight_type(V)}()
    for wv in basis(V; as_weight_vectors=true)
        push!(ws, wv)
    end
    return ws
end

function to_string(V::HighestWeightModule)
    return "V" * vec_subscript(highest_weight(V).weight)
end