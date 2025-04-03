export WeightSpace,
    space

struct WeightSpace{T<:AbstractSpace, W<:Weight}
    weight::W
    space::T
end

WeightSpace(weight::Vector, space::Vector) = WeightSpace(Weight(weight), MatrixVectorSpace(space))
WeightSpace{T,W}(wv::WeightVector) where {T,W} = WeightSpace{T,W}(weight(wv), T(vector(wv)))

Base.convert(
    ::Type{WeightSpace{T, W}},
    ws::WeightSpace
) where {T<:AbstractSpace, W<:Weight} = WeightSpace(
    convert(W, weight(ws)),
    convert(T, space(ws))
)

weight(ws::WeightSpace) = ws.weight
space(ws::WeightSpace) = ws.space
dim(ws::WeightSpace) = dim(space(ws))
field_type(::WeightSpace{T}) where {T} = field_type(T)

function Base.show(io::IO, ::MIME"text/plain", ws::WeightSpace)
    println(io, "WeightSpace of dimension $(dim(ws))")
    print(io, " weight: $(weight(ws))")
end

function Base.show(io::IO, ws::WeightSpace)
    print(io, "WeightSpace of dimension $(dim(ws))")
end

function Base.iterate(ws::WeightSpace, state=1)
    state > dim(ws) && return nothing
    return (WeightVector(weight(ws), basis(space(ws), state)), state+1)
end