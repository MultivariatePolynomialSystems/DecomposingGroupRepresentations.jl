export Weight,
    WeightVector,
    weight,
    vector,
    WeightSpace,
    space,
    WeightStructure,
    weight_space


struct Weight{T}
    weight::Vector{T}
end

Base.show(io::IO, ::MIME"text/plain", w::Weight) = print(io, "Weight: ", w.weight)
Base.show(io::IO, w::Weight) = print(io, w.weight)
Base.zero(::Type{Weight{T}}) where T = Weight(zero(T))
Base.zero(w::Weight) = Weight(zero(w.weight))
Base.vcat(w₁::Weight{T}, w₂::Weight{T}) where T = Weight(vcat(w₁.weight, w₂.weight))
Base.hash(w::Weight, h::UInt) = hash(w.weight, h)
Base.:(==)(w₁::Weight, w₂::Weight) = w₁.weight == w₂.weight

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

struct WeightSpace{F, T <: AbstractVectorSpace{F}, W<:Weight} # TODO: remove F?
    weight::W
    space::T
end

WeightSpace(weight::Vector, space::Vector) = WeightSpace(Weight(weight), MatrixVectorSpace(space))
WeightSpace{F,T,W}(wv::WeightVector) where {F,T,W} = WeightSpace{F,T,W}(weight(wv), T(vector(wv)))

Base.convert(
    ::Type{WeightSpace{F, T, W}},
    ws::WeightSpace
) where {F, T<:AbstractVectorSpace{F}, W<:Weight} = WeightSpace(
    convert(W, weight(ws)),
    convert(T, space(ws))
)

weight(ws::WeightSpace) = ws.weight
space(ws::WeightSpace) = ws.space
dim(ws::WeightSpace) = dim(space(ws))

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

struct WeightStructure{F, T<:AbstractVectorSpace{F}, W<:Weight}
    weights::Vector{W} # TODO: remove?
    dict::Dict{W, WeightSpace{F, T, W}}
end

WeightStructure{F,T,W}() where {F,T,W} = WeightStructure{F,T,W}(Weight[], Dict())

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

# TODO: improve
function WeightStructure(w_spaces::Vector{<:WeightSpace{F,T,W}}) where {F,T,W}
    ws = WeightStructure{F,T,W}()
    for w_space in w_spaces
        push!(ws, w_space)
    end
    return ws
end

Base.convert(
    ::Type{WeightStructure{F, T, W}},
    ws::WeightStructure
) where {F, T<:AbstractVectorSpace{F}, W<:Weight} = WeightStructure(
    convert(Vector{W}, ws.weights),
    convert(Dict{W, WeightSpace{F, T, W}}, ws.dict)
)

weights(ws::WeightStructure) = ws.weights
nweights(ws::WeightStructure) = length(ws.weights)
weight(ws::WeightStructure, i::Integer) = weights(ws)[i]
weight_space(ws::WeightStructure, i::Integer) = ws[weight(ws, i)]
dims(ws::WeightStructure) = [dim(space(ws[w])) for w in weights(ws)]
dim(ws::WeightStructure) = sum(dims(ws))
Base.length(ws::WeightStructure) = nweights(ws)
Base.getindex(ws::WeightStructure{F,T,W}, weight::W) where {F,T,W<:Weight} = ws.dict[weight]
Base.getindex(ws::WeightStructure, i::Integer) = ws[weight(ws, i)]
Base.setindex!(ws::WeightStructure, ws_new::WeightSpace, weight::Weight) = ws.dict[weight] = ws_new
Base.haskey(ws::WeightStructure, w::Weight) = haskey(ws.dict, w)
zero_weight(ws::WeightStructure) = zero(first(weights(ws)))
field_space(::WeightStructure{F, T}) where {F, T} = field_space(T)

function Base.show(io::IO, ::MIME"text/plain", ws::WeightStructure)
    println(io, "WeightStructure of $(dim(ws))-dimensional vector space")
    println(io, " $(nweights(ws)) weights: ", join(weights(ws), ", "))
    println(io, " dimensions of $(nweights(ws)) weight spaces: ", join(dims(ws), ", "))
    print(io, " max weight space dimension: ", maximum(dims(ws)))
end

function Base.show(io::IO, ws::WeightStructure)
    print(io, "WeightStructure of $(dim(ws))-dimensional vector space")
end

function Base.iterate(ws::WeightStructure, state=1)
    state > nweights(ws) && return nothing
    return (ws[state], state+1)
end

function Base.push!(ws::WeightStructure, w_space::WeightSpace)
    w = weight(w_space)
    if haskey(ws, w)
        ws[w] = WeightSpace(w, +(space(ws[w]), space(w_space)))
    else
        push!(weights(ws), w)
        ws[w] = w_space
    end
    return ws
end

Base.push!(ws::WeightStructure{F,T,W}, wv::WeightVector) where {F,T,W} = push!(ws, WeightSpace{F,T,W}(wv))

function sym_weight_structure(ws::WeightStructure{F}, d::Int) where F
    d == 0 && return WeightStructure([WeightSpace(zero_weight(ws), field_space(ws))])
    d == 1 && return ws
    combs = multiexponents(; degree=d, nvars=nweights(ws))
    new_weights_dict = Dict{Vector{Int}, Vector{typeof(combs[1])}}()
    for comb in combs
        w = sum([comb.nzval[i]*weight(ws, comb.nzind[i]) for i in 1:length(comb.nzind)])
        if haskey(new_weights_dict, w)
            push!(new_weights_dict[w], comb)
        else
            new_weights_dict[w] = [comb]
        end
    end
    new_weights = [zeros(Int, 0) for _ in 1:length(new_weights_dict)]
    new_weight_spaces = [zeros(ComplexF64, 0, 0) for _ in 1:length(new_weights_dict)]
    for (i, (weight, combs)) in enumerate(new_weights_dict)
        new_weights[i] = weight
        Ms = [⊙(weight_spaces(ws, comb.nzind; as_matrices=true), comb.nzval, mexps) for comb in combs]
        new_weight_spaces[i] = hcat(Ms...)
    end
    return WeightStructure(new_weights, new_weight_spaces)
end