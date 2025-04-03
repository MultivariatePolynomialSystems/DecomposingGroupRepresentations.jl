module DecomposingRepresentations

using SparseArrays: SparseMatrixCSC, sparse, findnz, spzeros
using LinearAlgebra: Diagonal, nullspace, I, dot, norm

import DynamicPolynomials
using DynamicPolynomials: @polyvar, Variable, Monomial, Polynomial, Commutative, CreationOrder, Graded, LexOrder, AbstractPolynomialLike, AbstractPolynomial
using DynamicPolynomials: subs, monomials, coefficients, differentiate, variables, nvariables, MonomialVector
export @polyvar, Variable, Monomial, Polynomial, Commutative, CreationOrder, Graded, LexOrder

using Combinatorics: partitions, multiset_permutations, combinations, with_replacement_combinations

using Base.Iterators: product, flatten
using Crayons

include("utils/basic.jl")
include("utils/Gauss-Jordan.jl")
include("types.jl")

include("vector-spaces/basic.jl")
include("vector-spaces/composite.jl")
include("vector-spaces/polynomials.jl")

include("lie-groups/lie-algebras/weights/weight.jl")
include("lie-groups/lie-algebras/weights/weight-vector.jl")
include("lie-groups/lie-algebras/weights/weight-space.jl")
include("lie-groups/lie-algebras/weights/weight-struct.jl")
include("lie-groups/lie-algebras/algebras.jl")
include("lie-groups/lie-algebras/elements.jl")
include("lie-groups/lie-algebras/weights/decompose-weights.jl")

include("lie-groups/lie-groups/groups.jl")
include("lie-groups/lie-groups/elements.jl")
include("lie-groups/lie-groups/actions/actions.jl")
include("lie-groups/lie-groups/actions/act.jl")
include("lie-groups/lie-groups/actions/action-ws.jl")

include("vector-spaces/hw-module.jl")
include("representations/irred-reprs.jl")
include("representations/isotypic-comps.jl")
include("representations/reprs.jl")

include("decompositions/irred-decomp.jl")
include("decompositions/isotypic-decomp.jl")

end
