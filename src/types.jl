export GroupType, Lie, Finite, Mixed
export AbstractGroup, AbstractDirectProductGroup, AbstractGroupElem, AbstractAction

abstract type GroupType end
struct Lie <: GroupType end
struct Finite <: GroupType end
struct Mixed <: GroupType end

abstract type AbstractGroup{T<:GroupType, F} end

lie_algebra(::AbstractGroup{Lie}) = error("Not implemented")


abstract type AbstractDirectProductGroup{T<:GroupType, F} <: AbstractGroup{T, F} end

abstract type AbstractGroupElem end

group(::AbstractGroupElem) = error("Not implemented")


abstract type AbstractAction{T<:GroupType, F} end

group(::AbstractAction) = error("Not implemented")


abstract type AbstractLieAlgebra{F} end

name(::AbstractLieAlgebra) = error("Not implemented")
basis(::AbstractLieAlgebra) = error("Not implemented")
chevalley_basis(::AbstractLieAlgebra) = error("Not implemented")
dim(::AbstractLieAlgebra) = error("Not implemented")
rank(::AbstractLieAlgebra) = error("Not implemented")


abstract type AbstractLieAlgebraElem end

algebra(::AbstractLieAlgebraElem) = error("Not implemented")


abstract type AbstractVectorSpace{F} end # TODO: remove F?

name(::AbstractVectorSpace) = error("Not implemented")
basis(::AbstractVectorSpace) = error("Not implemented")
dim(::AbstractVectorSpace) = error("Not implemented")


abstract type AbstractDirectSum{F} <: AbstractVectorSpace{F} end

spaces(::AbstractDirectSum) = error("Not implemented")


abstract type AbstractSymmetricPower{F} <: AbstractVectorSpace{F} end

base_space(::AbstractSymmetricPower) = error("Not implemented")
power(::AbstractSymmetricPower) = error("Not implemented")


abstract type AbstractTensorProduct{F} <: AbstractVectorSpace{F} end

spaces(::AbstractTensorProduct) = error("Not implemented")


abstract type AbstractRepresentation{T<:AbstractGroup, S<:AbstractVectorSpace} end

group(::AbstractRepresentation) = error("Not implemented")
action(::AbstractRepresentation) = error("Not implemented")
space(::AbstractRepresentation) = error("Not implemented")
degree(ρ::AbstractRepresentation) = dim(space(ρ))
irreducibles(::AbstractRepresentation) = error("Not implemented")
isotypic_components(::AbstractRepresentation) = error("Not implemented")
