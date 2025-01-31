export MatrixVectorSpace,
    VectorSpace


struct MatrixVectorSpace{T <: Number} <: AbstractVectorSpace
    basis::Matrix{T}
end


struct VectorSpace{T} <: AbstractVectorSpace
    basis::Vector{T}
end