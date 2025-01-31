abstract type AbstractReductiveGroup end

abstract type AbstractReductiveLieGroup <: AbstractReductiveGroup end

abstract type AbstractFiniteGroup <: AbstractReductiveGroup end

abstract type AbstractGroupAction end

group(::AbstractGroupAction) = error("Not implemented")


abstract type AbstractLieAlgebra end

name(::AbstractLieAlgebra) = error("Not implemented")
basis(::AbstractLieAlgebra) = error("Not implemented")
chevalley_basis(::AbstractLieAlgebra) = error("Not implemented")
dim(::AbstractLieAlgebra) = error("Not implemented")
rank(::AbstractLieAlgebra) = error("Not implemented")


abstract type AbstractLieAlgebraElem end

algebra(::AbstractLieAlgebraElem) = error("Not implemented")


abstract type AbstractVectorSpace end

name(::AbstractVectorSpace) = error("Not implemented")
basis(::AbstractVectorSpace) = error("Not implemented")
dim(::AbstractVectorSpace) = error("Not implemented")


abstract type AbstractDirectSum <: AbstractVectorSpace end

spaces(::AbstractDirectSum) = error("Not implemented")


abstract type AbstractSymmetricPower <: AbstractVectorSpace end

base_space(::AbstractSymmetricPower) = error("Not implemented")
power(::AbstractSymmetricPower) = error("Not implemented")


abstract type AbstractTensorProduct <: AbstractVectorSpace end

spaces(::AbstractTensorProduct) = error("Not implemented")


abstract type AbstractRepresentation{T<:AbstractReductiveGroup, S<:AbstractVectorSpace} end

group(::AbstractRepresentation) = error("Not implemented")
action(::AbstractRepresentation) = error("Not implemented")
space(::AbstractRepresentation) = error("Not implemented")
degree(ρ::AbstractRepresentation) = dim(space(ρ))
irreducibles(::AbstractRepresentation) = error("Not implemented")
isotypic_components(::AbstractRepresentation) = error("Not implemented")
