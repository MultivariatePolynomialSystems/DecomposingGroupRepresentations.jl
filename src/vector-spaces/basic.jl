export MatrixVectorSpace,
    VectorSpace,
    VariableSpace


struct MatrixVectorSpace{F} <: AbstractVectorSpace{F}
    matrix::Matrix{F}
end

MatrixVectorSpace(v::Vector) = MatrixVectorSpace(V2M(v))

matrix(V::MatrixVectorSpace) = V.matrix
basis(V::MatrixVectorSpace) = matrix(V) # TODO: change to eachcol?
basis(V::MatrixVectorSpace, i::Integer) = matrix(V)[:, i]
dim(V::MatrixVectorSpace) = size(basis(V), 2)
Base.iszero(V::MatrixVectorSpace) = dim(V) == 0

Base.convert(
    ::Type{MatrixVectorSpace{T}},
    V::MatrixVectorSpace
) where {T} = MatrixVectorSpace(convert(Matrix{T}, matrix(V)))

# TODO: Use generic nullspace
function Base.:∩(V₁::MatrixVectorSpace, V₂::MatrixVectorSpace)
    N = nullspace(hcat(matrix(V₁), matrix(V₂)))
    return MatrixVectorSpace(matrix(V₁)*N[1:dim(V₁), :]) # TODO: pick V with smaller dim
end

# TODO: add gens?
struct VectorSpace{T, F} <: AbstractVectorSpace{F}
    basis::Vector{T}
end

VectorSpace{T, F}() where {T, F} = VectorSpace{T, F}(Vector{T}())
VectorSpace{T, F}(v::T) where {T, F} = VectorSpace{T, F}([v])
VectorSpace(F::DataType, vs::Vector{T}) where {T} = VectorSpace{T, F}(vs)

basis(V::VectorSpace) = V.basis
basis(V::VectorSpace, i::Integer) = basis(V)[i]
dim(V::VectorSpace) = length(basis(V))
Base.iszero(V::VectorSpace) = dim(V) == 0
field_space(::Type{VectorSpace{T, F}}) where {T, F} = VectorSpace{T, F}(T(one(F)))

Base.push!(V::VectorSpace{T}, v::T) where T = push!(V.basis, v)

function Base.show(io::IO, V::VectorSpace{T, F}; indent::Int=0) where {T, F}
    println(io, " "^indent, "VectorSpace of dimension $(dim(V))")
    println(io, " "^indent, " element type: $(T)")
    println(io, " "^indent, " number type (or field): $(F)")
    print(io, " "^indent, " basis: ", join(map(repr, basis(V)), ", "))
end

Base.:+(
    Vs::VectorSpace{T, F}...
) where {T<:Variable, F} = VectorSpace{T,F}(collect(Set(vcat([basis(V) for V in Vs]...))))
Base.rand(V::VectorSpace{T, F}) where {T, F} = sum(rand(F, dim(V)) .* basis(V))

variables(V::VectorSpace{<:Variable}) = basis(V)

function Base.:∩(
    V₁::VectorSpace{<:AbstractPolynomialLike, F},
    V₂::VectorSpace{<:AbstractPolynomialLike, F}
) where F
    all_mons = monomials(basis(V₁), basis(V₂))
    M₁ = coeffs_matrix(basis(V₁), all_mons)
    M₂ = coeffs_matrix(basis(V₂), all_mons)
    N = nullspace(hcat(M₁, M₂))
    N = hcat([div_by_lowest_magnitude(N[:,i], 1e-8) for i in 1:size(N, 2)]...)
    sparsify!(N, 1e-8)
    Vᵢ = M₁*N[1:dim(V₁), :]
    return VectorSpace(F, [sum(c .* all_mons) for c in eachcol(Vᵢ)])
end

function zero_combinations(F::Vector{<:AbstractPolynomial})
    mons = monomials(F)
    M = coeffs_matrix(F, mons)
    return eachcol(nullspace(M))
end