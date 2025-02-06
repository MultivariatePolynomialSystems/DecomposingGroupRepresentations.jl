module DecomposingRepresentations

using SparseArrays: SparseMatrixCSC

using DynamicPolynomials: @polyvar
export @polyvar

include("utils/basic.jl")
include("types.jl")

include("reprs/spaces.jl")

include("lie/lie-algebras/weights.jl")
include("lie/lie-algebras/algebras.jl")
include("lie/lie-algebras/elements.jl")

include("lie/lie-groups/groups.jl")

end
