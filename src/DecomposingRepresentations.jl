module DecomposingRepresentations

using SparseArrays: SparseMatrixCSC, sparse, findnz, spzeros
using LinearAlgebra: Diagonal, nullspace, I, dot, norm

using SymEngine: Basic, free_symbols
const Expression = Basic

using Combinatorics: partitions, multiset_permutations, combinations, with_replacement_combinations

using Base.Iterators: product

include("utils/basic.jl")
include("utils/Gauss-Jordan.jl")
include("types.jl")

include("vector-spaces/basic.jl")
include("vector-spaces/composite.jl")

include("lie-groups/lie-algebras/weights.jl")
include("lie-groups/lie-algebras/algebras.jl")
include("lie-groups/lie-algebras/elements.jl")

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
