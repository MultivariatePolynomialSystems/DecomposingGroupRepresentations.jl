export Weight,
    WeightVector,
    WeightSpace,
    WeightStructure


struct Weight{T}
    weight::Vector{T}
end

# Base.convert(::Type{Weight{T}}, v::Vector{T}) where T = Weight(v)

struct WeightVector{W, T}
    weight::Weight{W} # TODO: change to W
    vector::T
end

weight(wv::WeightVector) = wv.weight
vector(wv::WeightVector) = wv.vector

struct WeightSpace{F, T <: AbstractVectorSpace{F}, W} # TODO: <: AbstractVectorSpace?
    weight::Weight{W}
    space::T
end

WeightSpace(weight::Vector, space::Vector) = WeightSpace(Weight(weight), MatrixVectorSpace(space))

Base.convert(
    ::Type{WeightSpace{F, T, W}},
    ws::WeightSpace
) where {F, T<:AbstractVectorSpace{F}, W} = WeightSpace(
    convert(Weight{W}, weight(ws)),
    convert(T, space(ws))
)

weight(ws::WeightSpace) = ws.weight
space(ws::WeightSpace) = ws.space
dim(ws::WeightSpace) = dim(space(ws))

struct WeightStructure{F, T<:AbstractVectorSpace{F}, W}
    weights::Vector{Weight{W}}
    dict::Dict{Weight{W}, WeightSpace{F, T, W}}
end

WeightStructure(
    weights::Vector{<:Weight},
    weight_spaces::Vector{<:MatrixVectorSpace}
) = WeightStructure(
    weights,
    Dict(zip(weights, [WeightSpace(w, V) for (w, V) in zip(weights, weight_spaces)]))
)

WeightStructure(
    weights::Vector{<:Vector},
    weight_vectors::Vector{<:Vector}
) = WeightStructure(
    [Weight(w) for w in weights],
    [MatrixVectorSpace(v) for v in weight_vectors]
)

Base.convert(
    ::Type{WeightStructure{F, T, W}},
    ws::WeightStructure
) where {F, T<:AbstractVectorSpace{F}, W} = WeightStructure(
    convert(Vector{Weight{W}}, ws.weights),
    convert(Dict{Weight{W}, WeightSpace{F, T, W}}, ws.dict)
)

weights(ws::WeightStructure) = ws.weights
nweights(ws::WeightStructure) = length(ws.weights)
Base.length(ws::WeightStructure) = nweights(ws)
Base.haskey(ws::WeightStructure, w::Weight) = haskey(ws.dict, w)