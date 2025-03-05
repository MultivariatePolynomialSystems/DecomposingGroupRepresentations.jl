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


struct VariableSpace{F} <: AbstractVectorSpace{F}
    variables::Vector{Expression}

    function VariableSpace(vars::AbstractVector{Expression})
        v = Vector{Expression}(undef, length(vars))
        for (i, var) in enumerate(vars)
            fs = free_symbols(var)
            @assert length(fs) == 1 && fs[1] == var
            v[i] = var
        end
        return new{F}(v)
    end
end

variables(V::VariableSpace) = V.variables
nvariables(V::VariableSpace) = length(variables(V))
basis(V::VariableSpace) = variables(V)
basis(V::VariableSpace, i::Integer) = basis(V)[i]
dim(V::VariableSpace) = length(basis(V))
Base.iszero(V::VariableSpace) = dim(V) == 0

# TODO: add gens?
struct ExpressionSpace{F} <: AbstractVectorSpace{F}
    gens::Vector{Expression}
end

ExpressionSpace{F}() where F = ExpressionSpace{F}(Vector{Expression}())
ExpressionSpace{F}(v::Expression) where F = ExpressionSpace{F}([v])

gens(V::ExpressionSpace) = V.gens
basis(V::ExpressionSpace) = gens(V) # FIXME
basis(V::ExpressionSpace, i::Integer) = basis(V)[i]
dim(V::ExpressionSpace) = length(basis(V))
Base.iszero(V::ExpressionSpace) = dim(V) == 0

function Base.show(io::IO, V::ExpressionSpace{F}; indent::Int=0) where F
    println(io, " "^indent, "ExpressionSpace of with $(dim(V)) generators")
    println(io, " "^indent, " number type (or field): $(F)")
    print(io, " "^indent, " generators: ", join(map(repr, basis(V)), ", "))
end

field_space(::Type{ExpressionSpace{F}}) where F = ExpressionSpace{F}(Expression(1))
variables(V::ExpressionSpace) = free_symbols(basis(V))
nvariables(V::ExpressionSpace) = length(variables(V))
Base.push!(V::ExpressionSpace, v::Expression) = push!(V.basis, v)

Base.rand(V::AbstractVectorSpace{F}) where F = sum(rand(F, dim(V)) .* basis(V))

# Base.convert(::Type{VectorSpace{T₁, F}}, V::VectorSpace{T₂, F}) where {T₁, T₂, F} = VectorSpace{T₁, F}(convert(Vector{T₁}, basis(V)))

Base.:+(
    Vs::ExpressionSpace{F}...
) where F = ExpressionSpace{F}(vcat([gens(V) for V in Vs]...))

Base.:*(
    Vs::Vector{ExpressionSpace{F}}
) where F = ExpressionSpace{F}([prod(fs) for fs in product([basis(V) for V in Vs]...)][:])
Base.:*(
    Vs::Vector{ExpressionSpace{F}},
    muls::Vector{Int}
) where F = *(vcat([fill(V, mul) for (V, mul) in zip(Vs, muls)]...))

function Base.:∩(
    V₁::ExpressionSpace{F},
    V₂::ExpressionSpace{F};
    tol::Real=1e-5
) where F
    
end
