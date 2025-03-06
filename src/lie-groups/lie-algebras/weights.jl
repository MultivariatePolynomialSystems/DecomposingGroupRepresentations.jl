export Weight,
    WeightVector,
    weight,
    vector,
    WeightSpace,
    space,
    WeightStructure,
    weight_space,
    sym


struct Weight{T}
    weight::Vector{T}
end

Base.show(io::IO, ::MIME"text/plain", w::Weight) = print(io, "Weight: ", w.weight)
Base.show(io::IO, w::Weight) = print(io, w.weight)
Base.zero(::Type{Weight{T}}, n) where T = Weight(zeros(T, n))
Base.zero(w::Weight) = Weight(zero(w.weight))
Base.vcat(w₁::Weight{T}, w₂::Weight{T}) where T = Weight(vcat(w₁.weight, w₂.weight))
Base.hash(w::Weight, h::UInt) = hash(w.weight, h)
Base.:(==)(w₁::Weight, w₂::Weight) = w₁.weight == w₂.weight
Base.:*(n::Number, w::Weight) = Weight(n*w.weight)
Base.:+(w₁::Weight{T}, w₂::Weight{T}) where {T} = Weight{T}(w₁.weight + w₂.weight)
Base.promote_rule(::Type{Weight{T}}, ::Type{Weight{S}}) where {T,S} = Weight{promote_type(T, S)}

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

struct WeightSpace{T <: AbstractVectorSpace, W<:Weight}
    weight::W
    space::T
end

WeightSpace(weight::Vector, space::Vector) = WeightSpace(Weight(weight), MatrixVectorSpace(space))
WeightSpace{T,W}(wv::WeightVector) where {T,W} = WeightSpace{T,W}(weight(wv), T(vector(wv)))

Base.convert(
    ::Type{WeightSpace{T, W}},
    ws::WeightSpace
) where {T<:AbstractVectorSpace, W<:Weight} = WeightSpace(
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

struct WeightStructure{T<:AbstractVectorSpace, W<:Weight}
    weights::Vector{W} # Ordering needed for sym_weight_structure
    dict::Dict{W, WeightSpace{T, W}} # TODO: new type for WeightSpace{F,T,W}? Do all spaces have to be of the same type?
end

WeightStructure{T,W}() where {T,W} = WeightStructure{T,W}(Weight[], Dict())

WeightStructure(
    weights::Vector{<:Weight},
    weight_spaces::Vector{<:AbstractVectorSpace}
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
function WeightStructure(w_spaces::Vector{<:WeightSpace{T,W}}) where {T,W}
    ws = WeightStructure{T,W}()
    for w_space in w_spaces
        push!(ws, w_space)
    end
    return ws
end

Base.convert(
    ::Type{WeightStructure{T, W}},
    ws::WeightStructure
) where {T<:AbstractVectorSpace, W<:Weight} = WeightStructure(
    convert(Vector{W}, ws.weights),
    convert(Dict{W, WeightSpace{T, W}}, ws.dict)
)

weights(ws::WeightStructure) = ws.weights
weights(ws::WeightStructure, inds...) = getindex(weights(ws), inds...)
nweights(ws::WeightStructure) = length(ws.weights)
weight(ws::WeightStructure, i::Integer) = weights(ws)[i] # TODO: remove?
weight_space(
    ws::WeightStructure,
    i::Integer;
    as_space::Bool=false
) = as_space ? space(ws[weight(ws, i)]) : ws[weight(ws, i)]
weight_spaces(
    ws::WeightStructure;
    as_spaces::Bool=false
) = as_spaces ? [space(ws[w]) for w in weights(ws)] : [ws[w] for w in weights(ws)]
weight_spaces(
    ws::WeightStructure,
    inds...;
    as_spaces::Bool=false
) = as_spaces ? [space(ws[w]) for w in weights(ws, inds...)] : [ws[w] for w in weights(ws, inds...)]
dims(ws::WeightStructure) = [dim(space(ws[w])) for w in weights(ws)]
dim(ws::WeightStructure) = sum(dims(ws))
Base.length(ws::WeightStructure) = nweights(ws)
Base.isempty(ws::WeightStructure) = isempty(weights(ws))
Base.getindex(ws::WeightStructure{T,W}, weight::W) where {T,W<:Weight} = ws.dict[weight]
Base.getindex(ws::WeightStructure, i::Integer) = ws[weight(ws, i)]
Base.setindex!(ws::WeightStructure, ws_new::WeightSpace, weight::Weight) = ws.dict[weight] = ws_new
Base.haskey(ws::WeightStructure, w::Weight) = haskey(ws.dict, w)
zero_weight(ws::WeightStructure) = zero(first(weights(ws)))
field_space(::WeightStructure{T}) where {T} = field_space(T)
field_type(::WeightStructure{T}) where {T} = field_type(T)

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

Base.push!(ws::WeightStructure{T,W}, wv::WeightVector) where {T,W} = push!(ws, WeightSpace{T,W}(wv))

function sym(ws::WeightStructure{T,W}, d::Int) where {T, W}
    d == 0 && return WeightStructure([WeightSpace(zero_weight(ws), field_space(ws))])
    d == 1 && return ws
    combs = multiexponents(; degree=d, nvars=nweights(ws))
    new_weights_dict = Dict{W, Vector{typeof(first(combs))}}()
    for comb in combs
        w = sum([comb.nzval[i]*weight(ws, comb.nzind[i]) for i in 1:length(comb.nzind)])
        if haskey(new_weights_dict, w)
            push!(new_weights_dict[w], comb)
        else
            new_weights_dict[w] = [comb]
        end
    end
    new_ws = WeightStructure{T,W}()
    for (weight, combs) in new_weights_dict
        w_sps = [*(weight_spaces(ws, comb.nzind; as_spaces=true), comb.nzval) for comb in combs]
        push!(new_ws, WeightSpace(weight, +(w_sps...)))
    end
    return new_ws
end

function tensor(ws₁::WeightStructure{T,W}, ws₂::WeightStructure{T,W}) where {T, W}
    combs = Base.Iterators.product(1:nweights(ws₁), 1:nweights(ws₂))
    new_weights_dict = Dict{W, Vector{typeof(first(combs))}}()
    for comb in combs
        w = weight(ws₁, comb[1]) + weight(ws₂, comb[2])
        val = get(new_weights_dict, w, nothing)
        if isnothing(val)
            new_weights_dict[w] = [comb]
        else
            push!(new_weights_dict[w], comb)
        end
    end
    new_ws = WeightStructure{T,W}()
    for (weight, combs) in new_weights_dict
        w_sps = [*([weight_space(ws₁, comb[1]; as_space=true), weight_space(ws₂, comb[2]; as_space=true)]) for comb in combs]
        total_ws = +(w_sps...)
        push!(new_ws, WeightSpace(weight, total_ws))
    end
    return new_ws
end