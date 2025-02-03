export MatrixVectorSpace,
    VectorSpace


struct MatrixVectorSpace{T} <: AbstractVectorSpace
    basis::Matrix{T}
end


struct VectorSpace{T} <: AbstractVectorSpace
    basis::Vector{T}
end