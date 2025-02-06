export LieAlgebraElem,
    RootElem


struct LieAlgebraElem{T<:AbstractReductiveLieAlgebra, S} <: AbstractLieAlgebraElem
    alg::T
    coeffs::Vector{S}
end

algebra(elem::LieAlgebraElem) = elem.alg
Base.size(elem::LieAlgebraElem) = size(algebra(elem))
Base.zero(alg::ReductiveLieAlgebra{T,W}) where {T,W} = LieAlgebraElem(alg, zeros(T, dim(alg)))
Base.zero(alg::ScalingLieAlgebra) = LieAlgebraElem(alg, zeros(Int, dim(alg)))
Base.randn(alg::ReductiveLieAlgebra{T,W}) where {T,W} = LieAlgebraElem(alg, randn(T, dim(alg)))
Base.randn(alg::ScalingLieAlgebra) = LieAlgebraElem(alg, randn(ComplexF64, dim(alg)))

# called by Shift+Enter
function Base.show(io::IO, mime::MIME"text/plain", elem::LieAlgebraElem)
    println(io, "LieAlgebraElem of $(name(algebra(elem))):")
    println(io, " matrix representation:")
    show(io, mime, as_matrix(elem)) # TODO: add offset
end

# called by print and inside vectors/matrices
function Base.show(io::IO, elem::LieAlgebraElem)
    print(io, "LieAlgebraElem of $(name(algebra(elem))) with coefficients: ")
    show(io, elem.coeffs)
end

Base.:*(a, elem::LieAlgebraElem) = LieAlgebraElem(elem.alg, a*elem.coeffs)
Base.:*(elem::LieAlgebraElem, a) = LieAlgebraElem(elem.alg, a*elem.coeffs)

function Base.:+(X::LieAlgebraElem, Y::LieAlgebraElem)
    @assert algebra(X) == algebra(Y)
    return LieAlgebraElem(algebra(X), X.coeffs+Y.coeffs)
end

function as_matrix(elem::LieAlgebraElem)
    return sum([elem.coeffs[i]*Bᵢ for (i, Bᵢ) in enumerate(basis(algebra(elem); as_matrices=true))])
end


struct SumLieAlgebraElem <: AbstractLieAlgebraElem
    alg::SumReductiveLieAlgebra
    elems::Vector{LieAlgebraElem}
end

algebra(elem::SumLieAlgebraElem) = elem.alg
elems(elem::SumLieAlgebraElem) = elem.elems
Base.randn(alg::SumReductiveLieAlgebra) = SumLieAlgebraElem(alg, [randn(a) for a in algebras(alg)])
Base.zero(alg::SumReductiveLieAlgebra) = SumLieAlgebraElem(alg, [zero(a) for a in algebras(alg)])
Base.getindex(elem::SumLieAlgebraElem, i::Int) = elems(elem)[i]
Base.setindex!(elem::SumLieAlgebraElem, val, i) = (elems(elem)[i] = val)


struct RootElem{T<:AbstractLieAlgebraElem} <: AbstractLieAlgebraElem
    elem::T
    root::Root
end

algebra(elem::RootElem) = algebra(elem.elem)
root(elem::RootElem) = elem.root

act(
    elem::RootElem,
    wv::WeightVector,
    a::AbstractGroupAction
) = WeightVector(root(elem)+weight(wv), act(elem.elem, vector(wv), a))
