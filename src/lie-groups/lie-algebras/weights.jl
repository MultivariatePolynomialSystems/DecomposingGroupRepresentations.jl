export Weight,
    WeightVector,
    WeightSpace,
    WeightStructure


struct Weight{T <: Number}
    weight::Vector{T}
end


struct WeightVector{T <: Number, S}
    weight::Weight{T}
    vector::S
end

weight(wv::WeightVector) = wv.weight
vector(wv::WeightVector) = wv.vector

struct WeightSpace{T <: Number, S <: AbstractVectorSpace} # TODO: <: AbstractVectorSpace?
    weight::Weight{T}
    space::S
end

weight(ws::WeightSpace) = ws.weight
space(ws::WeightSpace) = ws.space
dim(ws::WeightSpace) = dim(space(ws))

struct WeightStructure{T <: Number, S <: AbstractVectorSpace}
    weights::Vector{Weight{T}}
    dict::Dict{Weight{T}, WeightSpace{T, S}}
end

weights(ws::WeightStructure) = ws.weights
nweights(ws::WeightStructure) = length(ws.weights)
Base.length(ws::WeightStructure) = nweights(ws)
Base.haskey(ws::WeightStructure, w::Weight) = haskey(ws.dict, w)