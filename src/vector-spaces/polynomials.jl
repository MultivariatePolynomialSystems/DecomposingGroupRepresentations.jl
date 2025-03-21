export FixedDegreePolynomials,
    FixedMultidegreePolynomials

struct FixedDegreePolynomials{F,V,M} <: AbstractSymmetricPower{Polynomial{V,M,F}, F}
    variables::Vector{Variable{V,M}}
    degree::Int
end

FixedDegreePolynomials(
    vars::Vector{Variable{V,M}},
    d::Int
) where {V,M} = FixedDegreePolynomials{ComplexF64,V,M}(unique(vars), d)

var_space(V::FixedDegreePolynomials{F}) where F = VariableSpace{F}(V.variables)
degree(V::FixedDegreePolynomials) = V.degree
dim(V::FixedDegreePolynomials) = num_mons(nvariables(V), degree(V))
DynamicPolynomials.variables(V::FixedDegreePolynomials) = V.variables
DynamicPolynomials.nvariables(V::FixedDegreePolynomials) = length(V.variables)

function Base.show(io::IO, V::FixedDegreePolynomials)
    println(io, "FixedDegreePolynomials space of dimension $(dim(V))")
    println(io, " variables: ", join(map(repr, variables(V)), ", "))
    print(io, " degree: ", V.degree)
end


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
dim(V::FixedMultidegreePolynomials) = prod([num_mons(length(vars), d) for (vars, d) in zip(var_groups(V), degrees(V))])
var_spaces(V::FixedMultidegreePolynomials{F}) where F = [VariableSpace{F}(vars) for vars in V.var_groups]
DynamicPolynomials.variables(V::FixedMultidegreePolynomials) = vcat(V.var_groups...)
DynamicPolynomials.nvariables(V::FixedMultidegreePolynomials) = sum(length(vars) for vars in V.var_groups)

function Base.show(io::IO, V::FixedMultidegreePolynomials)
    println(io, "FixedMultidegreePolynomials space of dimension $(dim(V))")
    println(io, " variable groups: ", join([join(map(repr, var_group), ", ") for var_group in var_groups(V)], ", "))
    print(io, " multidegree: ", V.degrees)
end
