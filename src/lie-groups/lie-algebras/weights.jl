export Weight,
    WeightVector,
    WeightSpace,
    WeightStructure


struct Weight{T}
    weight::Vector{T}
end

# Base.convert(::Type{Weight{T}}, v::Vector{T}) where T = Weight(v)

struct WeightVector{T, W}
    weight::Weight{W}
    vector::T
end

weight(wv::WeightVector) = wv.weight
vector(wv::WeightVector) = wv.vector

struct WeightSpace{T <: AbstractVectorSpace, W} # TODO: <: AbstractVectorSpace?
    weight::Weight{W}
    space::T
end

WeightSpace{T, W}(
    w::Vector{W},
    v::Vector
) where {T <: AbstractVectorSpace, W} = WeightSpace(Weight(w), T(V2M(v)))

WeightSpace(
    w::Vector{W},
    v::Vector
) where {W} = WeightSpace{MatrixVectorSpace, W}(w, v)

weight(ws::WeightSpace) = ws.weight
space(ws::WeightSpace) = ws.space
dim(ws::WeightSpace) = dim(space(ws))

struct WeightStructure{T<:AbstractVectorSpace, W}
    weights::Vector{Weight{W}}
    dict::Dict{Weight{W}, WeightSpace{T, W}}
end

WeightStructure{T, W}(
    weights::Vector{Weight{W}},
    weight_spaces::Vector{<:Matrix}
) where {T <: AbstractVectorSpace, W} = WeightStructure(
    weights,
    Dict(zip(weights, [WeightSpace(w, T(M)) for (w, M) in zip(weights, weight_spaces)]))
)

WeightStructure(
    weights::Vector{Vector{W}},
    weight_vectors::Vector{<:Vector}
) where {W} = WeightStructure{MatrixVectorSpace, W}([Weight(w) for w in weights], [V2M(v) for v in weight_vectors])

weights(ws::WeightStructure) = ws.weights
nweights(ws::WeightStructure) = length(ws.weights)
Base.length(ws::WeightStructure) = nweights(ws)
Base.haskey(ws::WeightStructure, w::Weight) = haskey(ws.dict, w)