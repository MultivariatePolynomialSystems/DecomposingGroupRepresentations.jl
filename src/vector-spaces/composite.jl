export DirectSum,
    SymmetricPower,
    base_space,
    power,
    TensorProduct,
    factors


struct DirectSum{T, F, S<:AbstractSpace{T, F}} <: AbstractDirectSum{T, F}
    summands::Vector{S}
end

summands(V::DirectSum) = V.summands
nsummands(V::DirectSum) = length(summands(V))
dim(V::DirectSum) = sum(dim.(summands(V)))
basis(V::DirectSum) = vcat([basis(Vi) for Vi in summands(V)]...)


struct SymmetricPower{T, F, S<:AbstractSpace{T, F}} <: AbstractSymmetricPower{T, F}
    base_space::S
    power::Int

    function SymmetricPower{T, F, S}(base_space::S, power::Int) where {T, F, S<:AbstractSpace{T, F}}
        power < 1 && throw(ArgumentError("Power must be a positive integer"))
        power == 1 && return base_space
        return new(base_space, power)
    end
end

SymmetricPower(
    base_space::S,
    power::Int
) where {T, F, S<:AbstractSpace{T, F}} = SymmetricPower{T, F, S}(base_space, power)

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

to_string(V::SymmetricPower) = "Sym$(superscript(power(V)))($(to_string(base_space(V))))"


struct TensorProduct{T, F, S₁<:AbstractSpace{T,F}, S₂<:AbstractSpace{T,F}} <: AbstractTensorProduct{T, F}
    V₁::S₁
    V₂::S₂

    function TensorProduct{T,F,S₁,S₂}(V₁::S₁, V₂::S₂) where {T, F, S₁<:AbstractSpace{T,F}, S₂<:AbstractSpace{T,F}}
        # all_vars = variables(V₁) ∪ variables(V₂)
        # if length(all_vars) != nvariables(V₁) + nvariables(V₂)
        #     println("V₁: ", variables(V₁))
        #     println("V₂: ", variables(V₂))
        #     throw(ArgumentError("Spaces must have disjoint variables"))
        # end
        return new(V₁, V₂)
    end
end

TensorProduct(
    V₁::S₁,
    V₂::S₂
) where {T, F, S₁<:AbstractSpace{T,F}, S₂<:AbstractSpace{T,F}} = TensorProduct{T, F, S₁, S₂}(V₁, V₂)

function TensorProduct(
    factors::Vector{S}
) where {T, F, S<:AbstractSpace{T, F}}
    length(factors) == 1 && return factors[1]
    return TensorProduct(factors[1], TensorProduct(factors[2:end]))
end

factors(V::TensorProduct) = [V.V₁, V.V₂]
Base.first(V::TensorProduct) = V.V₁
second(V::TensorProduct) = V.V₂
dim(V::TensorProduct) = dim(V.V₁) * dim(V.V₂)

function Base.show(io::IO, V::TensorProduct{T, F}; indent::Int=0) where {T,F}
    println(io, " "^indent, "TensorProduct V₁ ⊗ V₂ of dimension $(dim(V))")
    println(io, " "^indent, " element type: ", T)
    print(io, " "^indent, " number type (or field): ", F)
end

to_string(V::TensorProduct) = join([to_string(Vᵢ) for Vᵢ in factors(V)], " ⊗ ")