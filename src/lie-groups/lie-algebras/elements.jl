export LieAlgebraElem,
    RootElem


struct LieAlgebraElem{T <: Number, S <: Number} <: AbstractLieAlgebraElem
    alg::ReductiveLieAlgebra{T, S}
    coeffs::Vector{S}
end