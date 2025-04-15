export FixedDegreePolynomials,
    FixedMultidegreePolynomials

@doc raw"""
    FixedDegreePolynomials <: AbstractSymmetricPower

A type representing a space of polynomials of a fixed degree. I.e., for the variables ``x_1, \dots, x_n`` and degree ``d`` gives
``\mathbb{F}[x_1, \dots, x_n]_{d}``, the space of polynomials in ``n`` variables of degree ``d``.

# Constructors
```julia
FixedDegreePolynomials(vars::Vector{<:Variable}, degree::Int)
```

# Examples
```jldoctest
julia> @polyvar x[1:3];

julia> V = FixedDegreePolynomials(x, 2)
FixedDegreePolynomials space of dimension 6
 variables: x₁, x₂, x₃
 degree: 2

julia> base_space(V)
VariableSpace with 3 variables
 number type (or field): ComplexF64
 variables: x₁, x₂, x₃

julia> power(V)
2
```
"""
struct FixedDegreePolynomials{F,V,M} <: AbstractSymmetricPower{Polynomial{V,M,F}, F}
    variables::Vector{Variable{V,M}}
    degree::Int
end

FixedDegreePolynomials(
    vars::Vector{Variable{V,M}},
    d::Int
) where {V,M} = FixedDegreePolynomials{ComplexF64,V,M}(unique(vars), d)

var_space(V::FixedDegreePolynomials{F}) where F = VariableSpace{F}(V.variables)
base_space(V::FixedDegreePolynomials) = var_space(V)
degree(V::FixedDegreePolynomials) = V.degree
power(V::FixedDegreePolynomials) = degree(V)
DynamicPolynomials.variables(V::FixedDegreePolynomials) = V.variables
DynamicPolynomials.nvariables(V::FixedDegreePolynomials) = length(V.variables)

function Base.show(io::IO, ::MIME"text/plain", V::FixedDegreePolynomials)
    println(io, "FixedDegreePolynomials space of dimension $(dim(V))")
    println(io, " variables: ", join(map(repr, variables(V)), ", "))
    print(io, " degree: ", V.degree)
end

function Base.show(io::IO, V::FixedDegreePolynomials)
    print(io, "FixedDegreePolynomials space of dimension $(dim(V))")
end


@doc raw"""
    FixedMultidegreePolynomials <: AbstractTensorProduct

A type representing a space of polynomials of a fixed multidegree in certain groups of variables. For example,
for 2 groups of variables ``x_1,\dots,x_m`` and ``y_1,\dots,y_n``, and degrees ``d_1`` and ``d_2``, this gives
```math
\mathbb{F}[x_1, \dots, x_m, y_1, \dots, y_n]_{(d_1, d_2)} \cong \mathbb{F}[x_1, \dots, x_m]_{d_1} \otimes \mathbb{F}[y_1, \dots, y_n]_{d_2},
```
the space of multi-homogeneous polynomials in ``m+n`` variables of degree ``d_1`` in the first group and 
``d_2`` in the second group.

# Constructors
```julia
FixedMultidegreePolynomials(var_groups::Vector{Vector{<:Variable}}, degrees::Vector{Int})
```

# Examples
```jldoctest
julia> @polyvar x[1:3] y[1:2];

julia> V = FixedMultidegreePolynomials([x, y], [1, 1])
FixedMultidegreePolynomials space of dimension 6
 variable groups: [x₁, x₂, x₃], [y₁, y₂]
 multidegree: [1, 1]

julia> fs = factors(V)
2-element Vector{FixedDegreePolynomials{ComplexF64, Commutative{CreationOrder}, Graded{LexOrder}}}:
 FixedDegreePolynomials space of dimension 3
 FixedDegreePolynomials space of dimension 2

julia> fs[2]
FixedDegreePolynomials space of dimension 2
 variables: y₁, y₂
 degree: 1
```
"""
struct FixedMultidegreePolynomials{F,V,M} <: AbstractTensorProduct{Polynomial{V,M,F}, F}
    var_groups::Vector{Vector{Variable{V,M}}}
    degrees::Vector{Int}
end

FixedMultidegreePolynomials(
    var_groups::Vector{Vector{Variable{V,M}}},
    degrees::Vector{Int}
) where {V,M} = FixedMultidegreePolynomials{ComplexF64,V,M}(var_groups, degrees)


var_groups(V::FixedMultidegreePolynomials) = V.var_groups
degrees(V::FixedMultidegreePolynomials) = V.degrees
factors(V::FixedMultidegreePolynomials) = [FixedDegreePolynomials(vars, d) for (vars, d) in zip(var_groups(V), degrees(V))]
dim(V::FixedMultidegreePolynomials) = prod([num_mons(length(vars), d) for (vars, d) in zip(var_groups(V), degrees(V))])
var_spaces(V::FixedMultidegreePolynomials{F}) where F = [VariableSpace{F}(vars) for vars in V.var_groups]
DynamicPolynomials.variables(V::FixedMultidegreePolynomials) = vcat(V.var_groups...)
DynamicPolynomials.nvariables(V::FixedMultidegreePolynomials) = sum(length(vars) for vars in V.var_groups)

function Base.show(io::IO, ::MIME"text/plain", V::FixedMultidegreePolynomials)
    println(io, "FixedMultidegreePolynomials space of dimension $(dim(V))")
    println(io, " variable groups: ", join(["[$(join(map(repr, var_group), ", "))]" for var_group in var_groups(V)], ", "))
    print(io, " multidegree: ", V.degrees)
end

function Base.show(io::IO, V::FixedMultidegreePolynomials)
    print(io, "FixedMultidegreePolynomials space of dimension $(dim(V))")
end