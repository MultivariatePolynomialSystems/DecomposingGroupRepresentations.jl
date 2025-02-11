module DecomposingRepresentations

using SparseArrays: SparseMatrixCSC, sparse, findnz
using LinearAlgebra: Diagonal

using DynamicPolynomials: @polyvar, Variable, Monomial, Polynomial, subs, Commutative, CreationOrder, Graded, LexOrder, AbstractPolynomialLike
export @polyvar, Variable, Monomial, Polynomials, Commutative, CreationOrder, Graded, LexOrder

include("utils/basic.jl")
include("types.jl")

include("reprs/spaces.jl")

include("lie/lie-algebras/weights.jl")
include("lie/lie-algebras/algebras.jl")
include("lie/lie-algebras/elements.jl")

include("lie/lie-groups/groups.jl")
include("lie/lie-groups/elements.jl")
include("lie/lie-groups/actions.jl")

end
