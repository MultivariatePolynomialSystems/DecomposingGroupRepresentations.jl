export LieAlgebraElem,
    RootElem,
    element,
    root


struct LieAlgebraElem{F, T<:AbstractLieAlgebra{F}} <: AbstractLieAlgebraElem
    alg::T
    coeffs::Vector{F}
end

function LieAlgebraElem(alg::T, coeffs::AbstractVector{F}) where {F, T<:AbstractLieAlgebra{F}}
    @assert length(coeffs) == dim(alg)
    return LieAlgebraElem{F, T}(alg, coeffs)
end

algebra(elem::LieAlgebraElem) = elem.alg
matrix(elem::LieAlgebraElem) = sum([elem.coeffs[i]*Bᵢ for (i, Bᵢ) in enumerate(basis(algebra(elem); as_matrices=true))])
Base.size(elem::LieAlgebraElem) = size(algebra(elem))
Base.zero(alg::LieAlgebra{F}) where F = LieAlgebraElem(alg, zeros(F, dim(alg)))
Base.zero(alg::ScalingLieAlgebra{F}) where F = LieAlgebraElem(alg, zeros(F, dim(alg)))
Base.iszero(elem::LieAlgebraElem) = iszero(elem.coeffs)
Base.rand(alg::LieAlgebra{F}) where F = LieAlgebraElem(alg, rand(F, dim(alg)))
Base.rand(alg::ScalingLieAlgebra{F}) where F = LieAlgebraElem(alg, rand(F, dim(alg)))

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


struct SumLieAlgebraElem{F, T<:LieAlgebraElem{F}} <: AbstractLieAlgebraElem
    alg::SumLieAlgebra{F}
    elems::Vector{T}
end

algebra(elem::SumLieAlgebraElem) = elem.alg
elements(elem::SumLieAlgebraElem) = elem.elems
Base.rand(alg::SumLieAlgebra) = SumLieAlgebraElem(alg, [rand(a) for a in algebras(alg)])
Base.zero(alg::SumLieAlgebra) = SumLieAlgebraElem(alg, [zero(a) for a in algebras(alg)])
Base.getindex(elem::SumLieAlgebraElem, i::Int) = elements(elem)[i]
Base.setindex!(elem::SumLieAlgebraElem, val, i) = (elements(elem)[i] = val)


struct RootElem{T<:AbstractLieAlgebraElem} <: AbstractLieAlgebraElem
    elem::T
    root::Root
end

RootElem(elem::AbstractLieAlgebraElem, root::Vector{Int}) = RootElem(elem, Root(root))

RootElem(
    alg::T,
    coeffs::Vector{F},
    root::Root
) where {F, T<:AbstractLieAlgebra{F}} = RootElem(LieAlgebraElem(alg, coeffs), root)

algebra(elem::RootElem) = algebra(elem.elem)
element(elem::RootElem) = elem.elem
root(elem::RootElem) = elem.root
matrix(elem::RootElem) = matrix(element(elem))
