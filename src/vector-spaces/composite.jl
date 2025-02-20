export DirectSum,
    SymmetricPowerSpace,
    SymmetricPowersSpace


struct DirectSumSpace{F, T<:AbstractVectorSpace{F}} <: AbstractDirectSum{F}
    spaces::Vector{T}
end

spaces(V::DirectSumSpace) = V.spaces
nspaces(V::DirectSumSpace) = length(spaces(V))
dim(V::DirectSumSpace) = sum(dim.(spaces(V)))


struct SymmetricPowerSpace{F, T<:AbstractVectorSpace{F}} <: AbstractSymmetricPower{F}
    V::T
    power::Int
end

dim(V::SymmetricPowerSpace) = num_mons(dim(V), power)


struct SymmetricPowersSpace{F, T<:AbstractVectorSpace{F}} <: AbstractDirectSum{F}
    space::T
    powers::Vector{Int}
end

SymmetricPowersSpace(
    V::AbstractVectorSpace{F},
    powers::AbstractVector{<:Integer}
) where F = SymmetricPowersSpace{F, typeof(V)}(V, collect(powers))

SymmetricPowersSpace(
    vars::Vector{<:Variable},
    powers::AbstractVector{<:Integer}
) = SymmetricPowersSpace(VariableSpace{ComplexF64}(vars), powers)

base_space(V::SymmetricPowersSpace) = V.space
powers(V::SymmetricPowersSpace) = V.powers
npowers(V::SymmetricPowersSpace) = length(powers(V))
dim(V::SymmetricPowersSpace) = sum([num_mons(dim(base_space(V)), p) for p in powers(V)])

function Base.show(io::IO, V::SymmetricPowersSpace; indent::Int=0)
    println(io, " "^indent, "SymmetricPowersSpace (direct sum of $(npowers(V)) symmetric powers)")
    println(io, " "^indent, " base space:")
    show(io, base_space(V); indent=2)
    println(io)
    print(io, " "^indent, " powers: $(powers(V))")
end