export Weight

struct Weight{T}
    weight::Vector{T}
end

Base.show(io::IO, ::MIME"text/plain", w::Weight) = print(io, "Weight: ", w.weight)
Base.show(io::IO, w::Weight) = print(io, w.weight)
Base.zero(::Type{Weight{T}}, n) where T = Weight(zeros(T, n))
Base.zero(w::Weight) = Weight(zero(w.weight))
Base.vcat(w₁::Weight{T}, w₂::Weight{T}) where T = Weight(vcat(w₁.weight, w₂.weight))
Base.vcat(ws::Weight{T}...) where T = Weight(vcat([w.weight for w in ws]...))
Base.hash(w::Weight, h::UInt) = hash(w.weight, h)
Base.:(==)(w₁::Weight, w₂::Weight) = w₁.weight == w₂.weight
Base.:*(n::Number, w::Weight) = Weight(n*w.weight)
Base.:+(w₁::Weight{T}, w₂::Weight{T}) where {T} = Weight{T}(w₁.weight + w₂.weight)
Base.promote_rule(::Type{Weight{T}}, ::Type{Weight{S}}) where {T,S} = Weight{promote_type(T, S)}
Base.getindex(w::Weight, i::Integer) = w.weight[i]
Base.getindex(w::Weight, inds...) = getindex(w.weight, inds...)
Base.length(w::Weight) = length(w.weight)
Base.eachindex(w::Weight) = eachindex(w.weight)

function str_irr(V::String, w::Weight)
    return V * vec_subscript(w.weight)
end

function str_irr_decomp(V::String, ws_with_muls::Dict{W, Int}) where W<:Weight
    return join([(mul==1 ? "" : "$(mul)") * str_irr(V, w) for (w, mul) in ws_with_muls], " ⊕ ")
end