export GroupType, Lie, Finite, Mixed

export AbstractGroup, AbstractDirectProductGroup, AbstractGroupElem, AbstractGroupAction
export algebra, group
export AbstractLieAlgebra, AbstractLieAlgebraElem
export name, basis, chevalley_basis, dim, rank
export AbstractVectorSpace, AbstractDirectSum, AbstractSymmetricPower, AbstractTensorProduct
export summands, nsummands, base_space, power, spaces
export AbstractGroupRepresentation
export action, space, irreducibles, isotypic_components


"""
    abstract type GroupType end

Abstract type representing a group type in the context of decomposing representations. The concrete types are [`Finite`](@ref), [`Lie`](@ref) and [`Mixed`](@ref).
"""
abstract type GroupType end


"""
    struct Lie <: GroupType end

Represents a Lie group type.
"""
struct Lie <: GroupType end


"""
    struct Finite <: GroupType end

Represents a finite group type. This type is used to categorize groups that have a finite number of elements.
"""
struct Finite <: GroupType end


"""
    struct Mixed <: GroupType end

Represents a mixed group type. This type is used in direct products of finite groups with Lie groups.
"""
struct Mixed <: GroupType end

"""
    AbstractGroup{T<:GroupType, F}

An abstract type representing a reductive group. The type `T` represents a [`GroupType`](@ref), while `F` represents the number field (or number type) over which the group is defined.
"""
abstract type AbstractGroup{T<:GroupType, F} end

"""
    algebra(::AbstractGroup{Lie, F}) -> AbstractLieAlgebra{F}

Returns the Lie algebra of a given Lie group.

# Examples
```julia-repl
julia> SO3 = LieGroup("SO", 3);

julia> algebra(SO3)
LieAlgebra 𝖘𝖔(3)
 number type (or field): ComplexF64
 weight type: Int64
 dimension: 3
 rank (dimension of Cartan subalgebra): 1
```
"""
algebra(::AbstractGroup{Lie}) = error("Not implemented")


"""
    AbstractDirectProductGroup{T<:GroupType, F} <: AbstractGroup{T, F}

An abstract type representing a direct product group.
"""
abstract type AbstractDirectProductGroup{T<:GroupType, F} <: AbstractGroup{T, F} end

# """
#     AbstractGroupElem

# An abstract type representing a group element.
# """
abstract type AbstractGroupElem end

# """
#     group(::AbstractGroupElem) -> AbstractGroup

# Returns the group to which the given group element belongs.
# """
group(::AbstractGroupElem) = error("Not implemented")


"""
    AbstractGroupAction{T<:GroupType, F}

An abstract type representing a group action. It is used to define how a group acts on a given vector space.
"""
abstract type AbstractGroupAction{T<:GroupType, F} end

"""
    group(::AbstractGroupAction) -> AbstractGroup

Returns the group associated with a given `AbstractGroupAction`.
"""
group(::AbstractGroupAction) = error("Not implemented")

"""
    algebra(a::AbstractGroupAction{Lie}) -> AbstractLieAlgebra

Returns the associated Lie algebra of the Lie group of the action.
"""
algebra(a::AbstractGroupAction{Lie}) = algebra(group(a))


"""
    AbstractLieAlgebra{F}

An abstract type representing a Lie algebra. The type `F` represents the number field (or number type) over which the Lie algebra is defined.
"""
abstract type AbstractLieAlgebra{F} end

"""
    name(::AbstractLieAlgebra) -> String

Returns the name of the given Lie algebra.
"""
name(::AbstractLieAlgebra) = error("Not implemented")

@doc raw"""
    basis(::AbstractLieAlgebra)

Returns a basis of the given Lie algebra. For example, the Lie algebra ``\mathfrak{so}(3, \mathbb{C})``,
consists of skew-symmetric matrices and hence its (standard) basis is given by the matrices
```math
\left\{ 
\begin{bmatrix} 0 & 0 & 0 \\ 0 & 0 & -1 \\ 0 & 1 & 0 \end{bmatrix},
\begin{bmatrix} 0 & 0 & 1 \\ 0 & 0 & 0 \\ -1 & 0 & 0 \end{bmatrix},
\begin{bmatrix} 0 & -1 & 0 \\ 1 & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}
\right\}
```
"""
basis(::AbstractLieAlgebra) = error("Not implemented")

@doc raw"""
    chevalley_basis(::AbstractLieAlgebra)

Returns the Chevalley basis of the given Lie algebra. For example, for the Lie algebra ``\mathfrak{so}(3, \mathbb{C})`` returns 
```math
\left\{ 
\underbrace{\begin{bmatrix} 0 & 0 & -1 \\ 0 & 0 & -i \\ 1 & i & 0 \end{bmatrix}}_{J_+},
\underbrace{\begin{bmatrix} 0 & 0 & 1 \\ 0 & 0 & -i \\ -1 & i & 0 \end{bmatrix}}_{J_-},
\underbrace{\begin{bmatrix} 0 & -i & 0 \\ i & 0 & 0 \\ 0 & 0 & 0 \end{bmatrix}}_{J_3}
\right\}
```
with the commutation relations
```math
\begin{aligned}
[J_3, J_+] &= J_+ \\
[J_3, J_-] &= -J_- \\
[J_+, J_-] &= 2J_3
\end{aligned}
```
"""
chevalley_basis(::AbstractLieAlgebra) = error("Not implemented")

"""
    dim(::AbstractLieAlgebra) -> Int

Returns the dimension of the given Lie algebra.
"""
dim(::AbstractLieAlgebra) = error("Not implemented")

"""
    rank(::AbstractLieAlgebra) -> Int

Returns the rank (i.e. the dimension of the Cartan subalgebra) of the given Lie algebra.
"""
rank(::AbstractLieAlgebra) = error("Not implemented")


# """
#     AbstractLieAlgebraElem

# An abstract type representing an element of a Lie algebra.
# """
abstract type AbstractLieAlgebraElem end

# """
#     algebra(::AbstractLieAlgebraElem) -> AbstractLieAlgebra

# Returns the Lie algebra to which the given Lie algebra element belongs.
# """
algebra(::AbstractLieAlgebraElem) = error("Not implemented")


"""
    AbstractVectorSpace{F}

An abstract type representing a vector space. The type `F` represents the number field (or number type) over which the vector space is defined.
"""
abstract type AbstractVectorSpace{T, F} end

field_type(::Type{<:AbstractVectorSpace{T, F}}) where {T, F} = F

"""
    basis(::AbstractVectorSpace)

Returns a basis of the given vector space.
"""
basis(::AbstractVectorSpace) = error("Not implemented")

"""
    dim(::AbstractVectorSpace) -> Int

Returns the dimension of the given vector space.
"""
dim(::AbstractVectorSpace) = error("Not implemented")


"""
    AbstractDirectSum{F} <: AbstractVectorSpace{F}

An abstract type representing a direct sum of vector spaces over the field (or number type) `F`.
"""
abstract type AbstractDirectSum{T, F} <: AbstractVectorSpace{T, F} end

"""
    summands(::AbstractDirectSum)

Returns the summands of the given direct sum of vector spaces.
"""
summands(::AbstractDirectSum) = error("Not implemented")

"""
    nsummands(::AbstractDirectSum) -> Int

Returns the number of summands in the given direct sum of vector spaces.
"""
nsummands(::AbstractDirectSum) = error("Not implemented")


"""
    AbstractSymmetricPower{F} <: AbstractVectorSpace{F}

An abstract type representing a symmetric power of a vector space over the field (or number type) `F`.
"""
abstract type AbstractSymmetricPower{T, F} <: AbstractVectorSpace{T, F} end

@doc raw"""
    base_space(::AbstractSymmetricPower) -> AbstractVectorSpace

Returns the base vector space of the given symmetric power. I.e. for ``\mathrm{Sym}^n(V)``, this function returns ``V``.
"""
base_space(::AbstractSymmetricPower) = error("Not implemented")

"""
    power(::AbstractSymmetricPower) -> Int

Returns the power of the given symmetric power.
"""
power(::AbstractSymmetricPower) = error("Not implemented")


"""
    AbstractTensorProduct{F} <: AbstractVectorSpace{F}

An abstract type representing a tensor product of vector spaces over the field `F`.
"""
abstract type AbstractTensorProduct{T, F} <: AbstractVectorSpace{T, F} end

"""
    spaces(::AbstractTensorProduct) -> Tuple{AbstractVectorSpace}

Returns the tuple of vector spaces that form the given tensor product.
"""
spaces(::AbstractTensorProduct) = error("Not implemented")


"""
    AbstractGroupRepresentation{T<:GroupType, S<:AbstractVectorSpace}

An abstract type representing a group representation. The type `T` represents a `GroupType`, and `S` represents an `AbstractVectorSpace`.
"""
abstract type AbstractGroupRepresentation{T<:GroupType, S<:AbstractVectorSpace} end

"""
    action(::AbstractGroupRepresentation) -> AbstractGroupAction

Returns the group action associated with the given group representation.
"""
action(::AbstractGroupRepresentation) = error("Not implemented")

"""
    group(ρ::AbstractGroupRepresentation) -> AbstractGroup

Returns the group associated with the given group representation.
"""
group(ρ::AbstractGroupRepresentation) = group(action(ρ))

"""
    space(::AbstractGroupRepresentation) -> AbstractVectorSpace

Returns the vector space on which the given group representation acts.
"""
space(::AbstractGroupRepresentation) = error("Not implemented")

"""
    dim(ρ::AbstractGroupRepresentation) -> Int

Returns the dimension of the vector space on which the given group representation acts.
"""
dim(ρ::AbstractGroupRepresentation) = dim(space(ρ))

"""
    irreducibles(::AbstractGroupRepresentation) -> Any

Returns the irreducible components of the given group representation.
"""
irreducibles(::AbstractGroupRepresentation) = error("Not implemented")

"""
    isotypic_components(::AbstractGroupRepresentation) -> Any

Returns the isotypic components of the given group representation.
"""
isotypic_components(::AbstractGroupRepresentation) = error("Not implemented")
