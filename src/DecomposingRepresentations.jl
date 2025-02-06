module DecomposingRepresentations

using SparseArrays: SparseMatrixCSC, sparse

using DynamicPolynomials: @polyvar, PolyVar
export @polyvar, PolyVar

include("utils/basic.jl")
include("types.jl")

include("reprs/spaces.jl")

include("lie/lie-algebras/weights.jl")
include("lie/lie-algebras/algebras.jl")
include("lie/lie-algebras/elements.jl")

include("lie/lie-groups/groups.jl")
include("lie/lie-groups/actions.jl")

end
