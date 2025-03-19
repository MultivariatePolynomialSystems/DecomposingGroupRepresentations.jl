export WeightVector,
    weight,
    vector

struct WeightVector{T, W<:Weight}
    weight::W
    vector::T
end

weight(wv::WeightVector) = wv.weight
vector(wv::WeightVector) = wv.vector

function Base.show(io::IO, ::MIME"text/plain", wv::WeightVector)
    println(io, "WeightVector with weight $(weight(wv))")
    print(io, " vector: $(repr(vector(wv)))")
end

function Base.show(io::IO, wv::WeightVector)
    print(io, "WeightVector with weight $(weight(wv))")
end