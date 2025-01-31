export ReductiveLieAlgebra


struct Root
    root::Vector{Int}
end

Base.:+(r::Root, w::Weight{T}) where {T <: Number} = Weight(r.root + w.weight)
Base.:+(w::Weight{T}, r::Root) where {T <: Number} = Weight(w.weight + r.root)

struct ChevalleyBasis{T <: Number}
    std_basis::Vector{Matrix{T}} # TODO: remove?
    cartan::Vector{Vector{T}} # given by coefficients in std_basis
    positive::Vector{Vector{T}} # given by coefficients in std_basis
    negative::Vector{Vector{T}} # given by coefficients in std_basis
    positive_roots::Vector{Root}
    negative_roots::Vector{Root}
end

struct ReductiveLieAlgebra{T <: Number, S <: Number} <: AbstractReductiveLieAlgebra
    name::String
    basis::ChevalleyBasis{S}
    weight_structure::WeightStructure{T, MatrixVectorSpace{S}}
    hw_spaces::Vector{WeightSpace{T, MatrixVectorSpace{S}}}
end

cartan_subalgebra(alg::ReductiveLieAlgebra) = [LieAlgebraElem(alg, coeffs) for coeffs in alg.basis.cartan]
positive_root_elements(
    alg::ReductiveLieAlgebra
) = [RootElem(alg, coeffs, root) for (coeffs, root) in zip(alg.basis.positive, alg.basis.positive_roots)]
negative_root_elements(
    alg::ReductiveLieAlgebra
) = [RootElem(alg, coeffs, root) for (coeffs, root) in zip(alg.basis.negative, alg.basis.negative_roots)]
weight_structure(alg::ReductiveLieAlgebra) = alg.weight_structure
weights(alg::ReductiveLieAlgebra) = weights(alg.weight_structure)
nweights(alg::ReductiveLieAlgebra) = nweights(alg.weight_structure)
hw_spaces(alg::ReductiveLieAlgebra) = alg.hw_spaces
