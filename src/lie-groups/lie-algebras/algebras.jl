export ReductiveLieAlgebra


struct ScalingLieAlgebra <: AbstractReductiveLieAlgebra
    name::String
    exps::SparseMatrixCSC{Int, Int} # every row is a vector u = [u₁,...,uₖ] which acts on vars by λᵘ
end

ScalingLieAlgebra(size::Int) = ScalingLieAlgebra("ℂ", sparse(ones(Int, 1, size)))

name(alg::ScalingLieAlgebra) = alg.name
dim(alg::ScalingLieAlgebra) = size(alg.exps, 1)
exponents(alg::ScalingLieAlgebra) = alg.exps
rank(alg::ScalingLieAlgebra) = dim(alg)
Base.size(alg::ScalingLieAlgebra) = size(alg.exps, 2)

function basis(alg::ScalingLieAlgebra; as_matrices::Bool=false)
    as_matrices && return [Diagonal(e) for e in eachrow(alg.exps)]
    coeffs = eye(ComplexF64, dim(alg))
    return [LieAlgebraElem(alg, c) for c in eachcol(coeffs)]
end

cartan_subalgebra(alg::ScalingLieAlgebra) = basis(alg)
positive_root_elements(::ScalingLieAlgebra) = LieAlgebraElem{ScalingLieAlgebra}[]
negative_root_elements(::ScalingLieAlgebra) = LieAlgebraElem{ScalingLieAlgebra}[]

struct Root
    root::Vector{Int}
end

Base.:+(r::Root, w::Weight{T}) where {T <: Number} = Weight(r.root + w.weight)
Base.:+(w::Weight{T}, r::Root) where {T <: Number} = Weight(w.weight + r.root)
Base.convert(::Type{Root}, v::Vector{Int}) = Root(v)

struct ChevalleyBasis{T}
    std_basis::Vector{Matrix{T}} # TODO: remove?
    cartan::Vector{Vector{T}} # given by coefficients in std_basis
    positive::Vector{Vector{T}} # given by coefficients in std_basis
    negative::Vector{Vector{T}} # given by coefficients in std_basis
    positive_roots::Vector{Root}
    negative_roots::Vector{Root}
end


struct ReductiveLieAlgebra{T, W} <: AbstractReductiveLieAlgebra
    name::String
    basis::ChevalleyBasis{T}
    weight_structure::WeightStructure{MatrixVectorSpace{T}, W}
    hw_spaces::Vector{WeightSpace{MatrixVectorSpace{T}, W}}
end

function so3_lie_algebra()
    X₁ = [0 0 0; 0 0 -1; 0 1 0]
    X₂ = [0 0 1; 0 0 0; -1 0 0]
    X₃ = [0 -1 0; 1 0 0; 0 0 0]
    cartan = [[0, 0, im]] # J₃ = 0*X₁ + 0*X₂ + im*X₃
    positive = [[im, -1, 0]] # J₊
    negative = [[im, 1, 0]] # J₋
    pos_roots = [[1]]
    neg_roots = [[-1]]
    ch_basis = ChevalleyBasis{Complex{Rational{Int}}}([X₁, X₂, X₃], cartan, positive, negative, pos_roots, neg_roots)
    ws = WeightStructure([[-1], [0], [1]], [[1, -im, 0], [0, 0, 1], [1, im, 0]])
    hw_spaces = [WeightSpace([1], [1, im, 0])]
    return ReductiveLieAlgebra("𝖘𝖔(3,ℂ)", ch_basis, ws, hw_spaces)
end

# TODO
function ReductiveLieAlgebra(name::String, size::Int)
    if name == "so" && size == 3
        return so3_lie_algebra()
    else
        error("Not implemented")
    end
end

name(alg::ReductiveLieAlgebra) = alg.name
dim(alg::ReductiveLieAlgebra) = length(alg.basis.std_basis)
rank(alg::ReductiveLieAlgebra) = length(alg.basis.cartan)
Base.size(alg::ReductiveLieAlgebra) = size(alg.basis.std_basis[1], 1)

function Base.show(io::IO, alg::ReductiveLieAlgebra)
    println(io, "ReductiveLieAlgebra $(name(alg))")
    println(io, " dimension: $(dim(alg))")
    print(io, " rank (dimension of Cartan subalgebra): $(rank(alg))")
end

function basis(alg::ReductiveLieAlgebra; as_matrices::Bool=false)
    if as_matrices
        return alg.basis.std_basis
    end
    coeffs = eye(ComplexF64, dim(alg))
    return [LieAlgebraElem(alg, c) for c in eachcol(coeffs)]
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


struct SumReductiveLieAlgebra <: AbstractReductiveLieAlgebra
    name::String
    algs::Vector{AbstractReductiveLieAlgebra}
end

name(alg::SumReductiveLieAlgebra) = alg.name
algebras(alg::SumReductiveLieAlgebra) = alg.algs
dim(alg::SumReductiveLieAlgebra) = sum([dim(a) for a in algebras(alg)])
rank(alg::SumReductiveLieAlgebra) = sum([rank(a) for a in algebras(alg)])

function Base.show(io::IO, alg::SumReductiveLieAlgebra)
    println(io, "SumReductiveLieAlgebra $(name(alg))")
    println(io, " dimension: $(dim(alg))")
    print(io, " rank (dimension of Cartan subalgebra): $(rank(alg))")
end

⊕(
    alg₁::AbstractReductiveLieAlgebra,
    alg₂::AbstractReductiveLieAlgebra
) = SumReductiveLieAlgebra("$(name(alg₁)) ⊕ $(name(alg₂))", [alg₁, alg₂])

⊕(
    alg₁::SumReductiveLieAlgebra,
    alg₂::AbstractReductiveLieAlgebra
) = SumReductiveLieAlgebra("$(name(alg₁)) ⊕ $(name(alg₂))", [algebras(alg₁)..., alg₂])

function get_elements(alg::SumReductiveLieAlgebra, sym::Symbol)
    elems = SumLieAlgebraElem[]
    for (i, a) in enumerate(algebras(alg))
        a_elems = eval(sym)(a)
        alg_elems = [zero(alg) for _ in a_elems]
        if sym == :positive_root_elements || sym == :negative_root_elements
            roots = [[zero(Int, rank(a)) for a in algebras(alg)] for _ in a_elems]
        end
        for (j, elem) in enumerate(a_elems)
            alg_elems[j][i] = elem
            if sym == :positive_root_elements || sym == :negative_root_elements
                roots[j][i] = root(elem)
                alg_elems[j] = RootElem(alg_elems[j], roots[j])
            end
        end
        append!(elems, alg_elems)
    end
    return elems
end

cartan_subalgebra(alg::SumReductiveLieAlgebra) = get_elements(alg, :cartan_subalgebra)
positive_root_elements(alg::SumReductiveLieAlgebra) = get_elements(alg, :positive_root_elements)
negative_root_elements(alg::SumReductiveLieAlgebra) = get_elements(alg, :negative_root_elements)
