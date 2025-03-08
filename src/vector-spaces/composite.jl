export DirectSum,
    SymmetricPower,
    base_space,
    power,
    SymmetricPowers,
    TensorProduct,
    factors,
    factor


struct DirectSum{F, S<:AbstractVectorSpace{F}} <: AbstractDirectSum{F}
    summands::Vector{S}
end

summands(V::DirectSum) = V.summands
nsummands(V::DirectSum) = length(summands(V))
dim(V::DirectSum) = sum(dim.(summands(V)))
basis(V::DirectSum) = vcat([basis(Vi) for Vi in summands(V)]...)


struct SymmetricPower{F, S<:AbstractVectorSpace{F}} <: AbstractSymmetricPower{F}
    base_space::S
    power::Int
end

SymmetricPower(
    vars::Vector{V},
    power::Integer
) where {V<:Variable} = SymmetricPower(VariableSpace{ComplexF64}(vars), power)

base_space(V::SymmetricPower) = V.base_space
power(V::SymmetricPower) = V.power
dim(V::SymmetricPower) = num_mons(dim(base_space(V)), power(V))
DynamicPolynomials.variables(V::SymmetricPower) = variables(base_space(V))
DynamicPolynomials.nvariables(V::SymmetricPower) = nvariables(base_space(V))

function Base.show(io::IO, V::SymmetricPower; indent::Int=0)
    println(io, " "^indent, "SymmetricPower of dimension $(dim(V))")
    println(io, " "^indent, " base space:")
    show(io, base_space(V); indent=indent+2)
    println(io)
    print(io, " "^indent, " power: $(power(V))")
end


struct SymmetricPowers{F, S<:AbstractVectorSpace{F}} <: AbstractDirectSum{F}
    space::S
    powers::Vector{Int}
end

SymmetricPowers(
    V::AbstractVectorSpace{F},
    powers::AbstractVector{<:Integer}
) where F = SymmetricPowers{F, typeof(V)}(V, collect(powers))

SymmetricPowers(
    vars::Vector{V},
    powers::AbstractVector{<:Integer}
) where {V<:Variable} = SymmetricPowers(VariableSpace{ComplexF64}(vars), powers)

base_space(V::SymmetricPowers) = V.space
powers(V::SymmetricPowers) = V.powers
dim(V::SymmetricPowers) = sum([num_mons(dim(base_space(V)), p) for p in powers(V)])
summands(V::SymmetricPowers) = [SymmetricPower(base_space(V), p) for p in powers(V)]
nsummands(V::SymmetricPowers) = length(powers(V))

function Base.show(io::IO, V::SymmetricPowers; indent::Int=0)
    println(io, " "^indent, "SymmetricPowersSpace (direct sum of $(nsummands(V)) symmetric powers)")
    println(io, " "^indent, " base space:")
    show(io, base_space(V); indent=2)
    println(io)
    print(io, " "^indent, " powers: ", join(powers(V), ", "))
end


struct TensorProduct{F, S<:AbstractVectorSpace{F}} <: AbstractTensorProduct{F}
    factors::Vector{S}

    function TensorProduct{F,S}(factors::Vector{S}) where {F, S<:AbstractVectorSpace{F}}
        all_vars = ∪([variables(V) for V in factors]...)
        sum_nvars = sum([nvariables(V) for V in factors])
        if length(all_vars) != sum_nvars
            throw(ArgumentError("Factors must have disjoint variables"))
        end
        return new(factors)
    end
end

TensorProduct(factors::Vector{S}) where {F, S<:AbstractVectorSpace{F}} = TensorProduct{F, S}(factors)

TensorProduct(
    V::AbstractVectorSpace{F},
    Vs::Vector{<:AbstractVectorSpace{F}}
) where F = TensorProduct(vcat([V], Vs))

factors(V::TensorProduct) = V.factors
factor(V::TensorProduct, i::Int) = factors(V)[i]
factors(V::TensorProduct, inds...) = getindex(factors(V), inds...)
nfactors(V::TensorProduct) = length(factors(V))
dim(V::TensorProduct) = prod([dim(Vᵢ) for Vᵢ in factors(V)])

function Base.show(io::IO, V::TensorProduct{F}; indent::Int=0) where F
    println(io, " "^indent, "TensorProduct of dimension $(dim(V)) ($(nfactors(V)) factors)")
    print(io, " "^indent, " number type (or field): ", F)
end