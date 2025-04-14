abstract type AbstractWeightStruct end

struct TensorWeightStruct <: AbstractWeightStruct
    w₁::Weight
    w₂::Weight
end

Base.hash(s::TensorWeightStruct, h::UInt) = Base.hash((s.w₁, s.w₂), h)
Base.:(==)(s₁::TensorWeightStruct, s₂::TensorWeightStruct) = ((s₁.w₁ == s₂.w₁) && (s₁.w₂ == s₂.w₂)) || ((s₁.w₁ == s₂.w₂) && (s₁.w₂ == s₂.w₁))

struct SymWeightStruct <: AbstractWeightStruct
    w::Weight
    p::Int
end

Base.hash(s::SymWeightStruct, h::UInt) = Base.hash((s.w, s.p), h)
Base.:(==)(s₁::SymWeightStruct, s₂::SymWeightStruct) = (s₁.w == s₂.w) && (s₁.p == s₂.p)
