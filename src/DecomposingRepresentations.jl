module DecomposingRepresentations

using SparseArrays: SparseMatrixCSC, sparse, findnz
using LinearAlgebra: Diagonal

using DynamicPolynomials: @polyvar, Variable, Monomial, Polynomial, subs, Commutative, CreationOrder, Graded, LexOrder, AbstractPolynomialLike
export @polyvar, Variable, Monomial, Polynomials, Commutative, CreationOrder, Graded, LexOrder

include("utils/basic.jl")
include("types.jl")

include("vector-spaces/basic.jl")
include("vector-spaces/composite.jl")

include("lie-groups/lie-algebras/weights.jl")
include("lie-groups/lie-algebras/algebras.jl")
include("lie-groups/lie-algebras/elements.jl")

include("lie-groups/lie-groups/groups.jl")
include("lie-groups/lie-groups/elements.jl")
include("lie-groups/lie-groups/actions.jl")

include("vector-spaces/hw-module.jl")
include("representations/irred-reprs.jl")
include("representations/isotypic-comps.jl")
include("representations/reprs.jl")

end
