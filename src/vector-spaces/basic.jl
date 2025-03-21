export MatrixVectorSpace,
    VariableSpace,
    PolySpace


struct MatrixVectorSpace{F} <: AbstractVectorSpace{Vector{F}, F}
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


struct VariableSpace{F, T<:Variable} <: AbstractVectorSpace{T, F}
    vars::Vector{T}

    VariableSpace{F, T}(vars::Vector{T}) where {F, T<:Variable} = new{F, T}(unique(vars))
end

VariableSpace{F}(var::T) where {F, T<:Variable} = VariableSpace{F, T}([var])
VariableSpace{F}(vars::Vector{T}) where {F, T<:Variable} = VariableSpace{F, T}(vars)

DynamicPolynomials.variables(V::VariableSpace) = V.vars
DynamicPolynomials.nvariables(V::VariableSpace) = length(variables(V))
basis(V::VariableSpace) = variables(V)
basis(V::VariableSpace, i::Integer) = basis(V)[i]
dim(V::VariableSpace) = nvariables(V)
Base.iszero(V::VariableSpace) = dim(V) == 0
Base.:+(V::VariableSpace{F}, W::VariableSpace{F}) where F = VariableSpace{F}(∪(variables(V), variables(W)))

function Base.show(io::IO, V::VariableSpace{F}; indent::Int=0) where F
    println(io, " "^indent, "VariableSpace with $(dim(V)) variables")
    println(io, " "^indent, " number type (or field): $(F)")
    print(io, " "^indent, " variables: ", join(map(repr, basis(V)), ", "))
end

Base.rand(V::AbstractVectorSpace{T, F}) where {T, F} = sum(rand(F, dim(V)) .* basis(V))


struct PolySpace{F, T<:Polynomial} <: AbstractVectorSpace{T, F}
    gens::Vector{T}
end

PolySpace{F, T}(v::T) where {F, T<:Polynomial} = PolySpace{F, T}([v])
PolySpace{F}(polys::Vector{T}) where {F, T<:Polynomial} = PolySpace{F, T}(polys)

gens(V::PolySpace) = V.gens
basis(V::PolySpace) = gens(V) # FIXME
basis(V::PolySpace, i::Integer) = basis(V)[i]
dim(V::PolySpace) = length(basis(V))
Base.iszero(V::PolySpace) = dim(V) == 0

function Base.show(io::IO, V::PolySpace{F}; indent::Int=0) where F
    println(io, " "^indent, "PolySpace with $(dim(V)) generators")
    println(io, " "^indent, " number type (or field): $(F)")
    print(io, " "^indent, " generators: ", join(map(repr, basis(V)), ", "))
end

field_space(V::PolySpace{F}) where F = PolySpace{F}(one(first(gens(V))))
DynamicPolynomials.variables(V::PolySpace) = ∪([variables(p) for p in gens(V)]...)
DynamicPolynomials.nvariables(V::PolySpace) = length(variables(V))
Base.push!(V::PolySpace{F, T}, v::T) where {F, T<:Polynomial} = push!(V.basis, v)
Base.convert(::Type{PolySpace{F}}, V::VariableSpace{F}) where F = PolySpace{F}(basis(V))

Base.:+(
    Vs::PolySpace{F}...
) where F = PolySpace{F}(vcat([gens(V) for V in Vs]...))

Base.:*(
    Vs::Vector{<:PolySpace{F}}
) where F = PolySpace{F}([prod(fs) for fs in product([basis(V) for V in Vs]...)][:])
Base.:*(
    Vs::Vector{<:PolySpace{F}},
    muls::Vector{Int}
) where F = *(vcat([fill(V, mul) for (V, mul) in zip(Vs, muls)]...))


function Base.:∩(
    V₁::AbstractVectorSpace{<:AbstractPolynomialLike, F},
    V₂::AbstractVectorSpace{<:AbstractPolynomialLike, F};
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
    return PolySpace{F}([sum(c .* all_mons) for c in eachcol(Vᵢ)])
end
