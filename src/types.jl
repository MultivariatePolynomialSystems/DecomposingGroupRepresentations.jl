abstract type AbstractReductiveGroup end

abstract type AbstractReductiveLieGroup <: AbstractReductiveGroup end

lie_algebra(::AbstractReductiveLieGroup) = error("Not implemented")


abstract type AbstractFiniteGroup <: AbstractReductiveGroup end

abstract type AbstractGroupElem end

group(::AbstractGroupElem) = error("Not implemented")


abstract type AbstractGroupAction end

group(::AbstractGroupAction) = error("Not implemented")

abstract type AbstractMatrixGroupAction end


abstract type AbstractReductiveLieAlgebra end # TODO: make parametric?

name(::AbstractReductiveLieAlgebra) = error("Not implemented")
basis(::AbstractReductiveLieAlgebra) = error("Not implemented")
chevalley_basis(::AbstractReductiveLieAlgebra) = error("Not implemented")
dim(::AbstractReductiveLieAlgebra) = error("Not implemented")
rank(::AbstractReductiveLieAlgebra) = error("Not implemented")


abstract type AbstractLieAlgebraElem end

algebra(::AbstractLieAlgebraElem) = error("Not implemented")


abstract type AbstractVectorSpace{F} end

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


abstract type AbstractRepresentation{F, T<:AbstractReductiveGroup, S<:AbstractVectorSpace{F}} end

group(::AbstractRepresentation) = error("Not implemented")
action(::AbstractRepresentation) = error("Not implemented")
space(::AbstractRepresentation) = error("Not implemented")
degree(ρ::AbstractRepresentation) = dim(space(ρ))
irreducibles(::AbstractRepresentation) = error("Not implemented")
isotypic_components(::AbstractRepresentation) = error("Not implemented")
