export GroupType, Lie, Finite, Mixed
export algebra, group
export AbstractGroup, AbstractDirectProductGroup, AbstractGroupElem, AbstractGroupAction
export basis
export action, dim

abstract type GroupType end
struct Lie <: GroupType end
struct Finite <: GroupType end
struct Mixed <: GroupType end

abstract type AbstractGroup{T<:GroupType, F} end

algebra(::AbstractGroup{Lie}) = error("Not implemented")


abstract type AbstractDirectProductGroup{T<:GroupType, F} <: AbstractGroup{T, F} end

abstract type AbstractGroupElem end

group(::AbstractGroupElem) = error("Not implemented")


abstract type AbstractGroupAction{T<:GroupType, F} end

group(::AbstractGroupAction) = error("Not implemented")
algebra(a::AbstractGroupAction{Lie}) = algebra(group(a))


abstract type AbstractLieAlgebra{F} end

name(::AbstractLieAlgebra) = error("Not implemented")
basis(::AbstractLieAlgebra) = error("Not implemented")
chevalley_basis(::AbstractLieAlgebra) = error("Not implemented")
dim(::AbstractLieAlgebra) = error("Not implemented")
rank(::AbstractLieAlgebra) = error("Not implemented")


abstract type AbstractLieAlgebraElem end

algebra(::AbstractLieAlgebraElem) = error("Not implemented")


abstract type AbstractVectorSpace{F} end # TODO: remove F?

basis(::AbstractVectorSpace) = error("Not implemented")
dim(::AbstractVectorSpace) = error("Not implemented")


abstract type AbstractDirectSum{F} <: AbstractVectorSpace{F} end

summands(::AbstractDirectSum) = error("Not implemented")
nsummands(::AbstractDirectSum) = error("Not implemented")


abstract type AbstractSymmetricPower{F} <: AbstractVectorSpace{F} end

base_space(::AbstractSymmetricPower) = error("Not implemented")
power(::AbstractSymmetricPower) = error("Not implemented")


abstract type AbstractTensorProduct{F} <: AbstractVectorSpace{F} end

spaces(::AbstractTensorProduct) = error("Not implemented")


abstract type AbstractGroupRepresentation{T<:GroupType, S<:AbstractVectorSpace} end

action(::AbstractGroupRepresentation) = error("Not implemented")
group(ρ::AbstractGroupRepresentation) = group(action(ρ))
space(::AbstractGroupRepresentation) = error("Not implemented")
space_basis(ρ::AbstractGroupRepresentation) = basis(space(ρ))
dim(ρ::AbstractGroupRepresentation) = dim(space(ρ))
irreducibles(::AbstractGroupRepresentation) = error("Not implemented")
isotypic_components(::AbstractGroupRepresentation) = error("Not implemented")
