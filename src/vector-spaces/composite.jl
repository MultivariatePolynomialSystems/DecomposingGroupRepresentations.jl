export DirectSum,
    SymmetricPowerSpace,
    base_space,
    power,
    SymmetricPowersSpace


struct DirectSumSpace{F, T<:AbstractVectorSpace{F}} <: AbstractDirectSum{F}
    summands::Vector{T}
end

summands(V::DirectSumSpace) = V.summands
nsummands(V::DirectSumSpace) = length(summands(V))
dim(V::DirectSumSpace) = sum(dim.(summands(V)))
basis(V::DirectSumSpace) = vcat([basis(Vi) for Vi in summands(V)]...)


struct SymmetricPowerSpace{F, T<:AbstractVectorSpace{F}} <: AbstractSymmetricPower{F}
    base_space::T
    power::Int
end

SymmetricPowerSpace(
    vars::Vector{V},
    power::Integer
) where {V<:Variable} = SymmetricPowerSpace(VectorSpace{V, ComplexF64}(vars), power)

base_space(V::SymmetricPowerSpace) = V.base_space
power(V::SymmetricPowerSpace) = V.power
dim(V::SymmetricPowerSpace) = num_mons(dim(base_space(V)), power(V))

function Base.show(io::IO, V::SymmetricPowerSpace; indent::Int=0)
    println(io, " "^indent, "SymmetricPowerSpace")
    println(io, " "^indent, " base space:")
    show(io, base_space(V); indent=indent+2)
    println(io)
    print(io, " "^indent, " power: $(power(V))")
end


struct SymmetricPowersSpace{F, T<:AbstractVectorSpace{F}} <: AbstractDirectSum{F}
    space::T
    powers::Vector{Int}
end

SymmetricPowersSpace(
    V::AbstractVectorSpace{F},
    powers::AbstractVector{<:Integer}
) where F = SymmetricPowersSpace{F, typeof(V)}(V, collect(powers))

SymmetricPowersSpace(
    vars::Vector{V},
    powers::AbstractVector{<:Integer}
) where {V<:Variable} = SymmetricPowersSpace(VectorSpace{V, ComplexF64}(vars), powers)

base_space(V::SymmetricPowersSpace) = V.space
powers(V::SymmetricPowersSpace) = V.powers
dim(V::SymmetricPowersSpace) = sum([num_mons(dim(base_space(V)), p) for p in powers(V)])
summands(V::SymmetricPowersSpace) = [SymmetricPowerSpace(base_space(V), p) for p in powers(V)]
nsummands(V::SymmetricPowersSpace) = length(powers(V))

function Base.show(io::IO, V::SymmetricPowersSpace; indent::Int=0)
    println(io, " "^indent, "SymmetricPowersSpace (direct sum of $(nsummands(V)) symmetric powers)")
    println(io, " "^indent, " base space:")
    show(io, base_space(V); indent=2)
    println(io)
    print(io, " "^indent, " powers: ", join(powers(V), ", "))
end