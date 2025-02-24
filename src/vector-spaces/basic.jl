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
VectorSpace(F::DataType, vars::Vector{T}) where T<:Variable = VectorSpace{T, F}(unique(vars))
VectorSpace(F::DataType, vs::Vector{T}) where T = VectorSpace{T, F}(vs)
VectorSpace(F::DataType, polys::Vector{T}) where T<:Polynomial = VectorSpace{T,F}(in_rref(polys))

basis(V::VectorSpace) = V.basis
basis(V::VectorSpace, i::Integer) = basis(V)[i]
dim(V::VectorSpace) = length(basis(V))

function Base.show(io::IO, V::VectorSpace{T, F}; indent::Int=0) where {T, F}
    println(io, " "^indent, "VectorSpace of dimension $(dim(V))")
    println(io, " "^indent, " element type: $(T)")
    println(io, " "^indent, " number type (or field): $(F)")
    print(io, " "^indent, " basis: ", join(map(repr, basis(V)), ", "))
end

field_space(::Type{VectorSpace{T, F}}) where {T, F} = VectorSpace{T, F}(one(T))
variables(V::VectorSpace{<:Variable}) = basis(V)
Base.iszero(V::VectorSpace) = dim(V) == 0
Base.rand(V::VectorSpace{T, F}) where {T, F} = sum(rand(F, dim(V)) .* basis(V))
Base.push!(V::VectorSpace{T}, v::T) where T = push!(V.basis, v)
Base.convert(::Type{VectorSpace{T₁, F}}, V::VectorSpace{T₂, F}) where {T₁, T₂, F} = VectorSpace{T₁, F}(convert(Vector{T₁}, basis(V)))

Base.:+(
    Vs::VectorSpace{T, F}...
) where {T<:Variable, F} = VectorSpace{T,F}(∪([basis(V) for V in Vs]...))
Base.:+(
    Vs::VectorSpace{T, F}...
) where {T<:Polynomial, F} = VectorSpace{T,F}(in_rref(vcat([basis(V) for V in Vs]...)))

Base.:*(
    Vs::Vector{VectorSpace{T, F}}
) where {T<:Polynomial, F} = VectorSpace{T, F}(in_rref([prod(fs) for fs in product([basis(V) for V in Vs]...)][:]))
Base.:*(
    Vs::Vector{VectorSpace{T, F}},
    muls::Vector{Int}
) where {T<:Polynomial, F} = *(vcat([fill(V, mul) for (V, mul) in zip(Vs, muls)]...))

function Base.:∩(
    V₁::VectorSpace{<:AbstractPolynomialLike, F},
    V₂::VectorSpace{<:AbstractPolynomialLike, F};
    tol::Real=1e-5
) where F
    all_mons = monomials(basis(V₁), basis(V₂))
    M₁ = coeffs_matrix(basis(V₁), all_mons)
    M₂ = coeffs_matrix(basis(V₂), all_mons)
    N = nullspace(hcat(M₁, M₂); atol=tol)
    size(N, 2) == 0 && return nothing
    N = hcat([div_by_lowest_magnitude(N[:,i], 1e-8) for i in 1:size(N, 2)]...)
    sparsify!(N, tol)
    Vᵢ = M₁*N[1:dim(V₁), :]
    return VectorSpace(F, [sum(c .* all_mons) for c in eachcol(Vᵢ)])
end

function zero_combinations(F::Vector{<:AbstractPolynomial}; tol::Real=1e-5)
    mons = monomials(F)
    M = coeffs_matrix(F, mons)
    return eachcol(nullspace(M; atol=tol))
end