export MatrixVectorSpace,
    VectorSpace,
    VariableSpace


struct MatrixVectorSpace{F} <: AbstractVectorSpace{F}
    basis::Matrix{F}
end

MatrixVectorSpace(v::Vector) = MatrixVectorSpace(V2M(v))

basis(V::MatrixVectorSpace) = V.basis
Base.convert(
    ::Type{MatrixVectorSpace{T}},
    V::MatrixVectorSpace
) where {T} = MatrixVectorSpace(convert(Matrix{T}, basis(V)))


struct VariableSpace{F} <: AbstractVectorSpace{F}
    vars::Vector{Variable}
end

basis(V::VariableSpace) = V.vars
dim(V::VariableSpace) = length(basis(V))

function Base.show(io::IO, V::VariableSpace{F}; indent::Int=0) where F
    println(io, " "^indent, "VariableSpace with $(dim(V)) variables")
    println(io, " "^indent, " variables: ", join(map(repr, V.vars), ", "))
    print(io, " "^indent, " number type (or field): $(F)")
end

Base.:+(
    Vs::VariableSpace{F}...
) where F = VariableSpace{F}(collect(Set(vcat([basis(V) for V in Vs]...))))

Base.rand(V::VariableSpace{F}) where F = sum(rand(F, length(basis(V))) .* basis(V))