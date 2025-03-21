export WeightStructure,
    weight_space,
    sym

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

function sym_weights_dict(ws::Vector{W}, p::Int) where W<:Weight
    combs = multiexponents(; degree=p, nvars=length(ws))
    ws_dict = Dict{W, Vector{typeof(first(combs))}}()
    for comb in combs
        w = sum([val*ws[ind] for (val, ind) in zip(comb.nzval, comb.nzind)])
        if haskey(ws_dict, w)
            push!(ws_dict[w], comb)
        else
            ws_dict[w] = [comb]
        end
    end
    return ws_dict
end

function sym_weight_struct(hw::W, A::AbstractLieAlgebra, p::Int, ws::WeightStructure{T,W}) where {T, W<:Weight}
    dir_sum_hws = sym(hw, A, p)
    ws_dict = sym_weights_dict(weights(ws), p)
    new_ws = WeightStructure{T,W}()
    for hw in dir_sum_hws
        combs = ws_dict[hw]
        w_sps = [*(weight_spaces(ws, comb.nzind; as_spaces=true), comb.nzval) for comb in combs]
        push!(new_ws, WeightSpace(hw, +(w_sps...)))
    end
    return new_ws
end

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
    for (w, combs) in new_weights_dict
        w_sps = [*(weight_spaces(ws, comb.nzind; as_spaces=true), comb.nzval) for comb in combs]
        push!(new_ws, WeightSpace(w, +(w_sps...)))
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