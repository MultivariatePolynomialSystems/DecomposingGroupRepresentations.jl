export MatrixVectorSpace,
    VectorSpace


struct MatrixVectorSpace{T} <: AbstractVectorSpace{T}
    basis::Matrix{T}
end

MatrixVectorSpace(v::Vector) = MatrixVectorSpace(V2M(v))

basis(V::MatrixVectorSpace) = V.basis
Base.convert(
    ::Type{MatrixVectorSpace{T}},
    V::MatrixVectorSpace
) where {T} = MatrixVectorSpace(convert(Matrix{T}, basis(V)))

struct VectorSpace{T} <: AbstractVectorSpace{T}
    basis::Vector{T}
end