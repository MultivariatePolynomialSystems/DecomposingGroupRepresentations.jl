export LieAlgebra


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

function show_basis(io::IO, alg::ScalingLieAlgebra; offset::Int=0)
    for i in 1:dim(alg)
        print(io, " "^offset, "[", join(alg.exps[i, :], ", "), "]")
        i < dim(alg) && print(io, '\n')
    end
end

function Base.show(io::IO, alg::ScalingLieAlgebra; offset::Int=0)
    println(io, " "^offset, "ScalingLieAlgebra $(name(alg))")
    println(io, " "^offset, " dimension: $(dim(alg))")
    println(io, " "^offset, " basis (diagonal matrices):")
    show_basis(io, alg, offset=offset+2)
end

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

function ChevalleyBasis(
    std_basis::Vector{Matrix{T1}},
    cartan::Vector{Vector{T2}},
    positive::Vector{Vector{T3}},
    negative::Vector{Vector{T4}},
    positive_roots::Vector{Vector{Int}},
    negative_roots::Vector{Vector{Int}}
) where {T1,T2,T3,T4}
    T = promote_type(T1, T2, T3, T4)
    return ChevalleyBasis{T}(std_basis, cartan, positive, negative, [Root(r) for r in positive_roots], [Root(r) for r in negative_roots])
end

Base.convert(
    ::Type{ChevalleyBasis{T}},
    ch_basis::ChevalleyBasis
) where {T} = ChevalleyBasis(
    convert(Vector{Matrix{T}}, ch_basis.std_basis),
    convert(Vector{Vector{T}}, ch_basis.cartan),
    convert(Vector{Vector{T}}, ch_basis.positive),
    convert(Vector{Vector{T}}, ch_basis.negative),
    ch_basis.positive_roots,
    ch_basis.negative_roots
)


struct LieAlgebra{F, W} <: AbstractReductiveLieAlgebra
    name::String
    basis::ChevalleyBasis{F}
    weight_structure::WeightStructure{F, MatrixVectorSpace{F}, W}
    hw_spaces::Vector{WeightSpace{F, MatrixVectorSpace{F}, W}}
end

function so3()
    X₁ = [0 0 0; 0 0 -1; 0 1 0]
    X₂ = [0 0 1; 0 0 0; -1 0 0]
    X₃ = [0 -1 0; 1 0 0; 0 0 0]
    cartan = [[0, 0, im]] # J₃ = 0*X₁ + 0*X₂ + im*X₃
    positive = [[im, -1, 0]] # J₊
    negative = [[im, 1, 0]] # J₋
    pos_roots = [[1]]
    neg_roots = [[-1]]
    ch_basis = ChevalleyBasis([X₁, X₂, X₃], cartan, positive, negative, pos_roots, neg_roots)
    ws = WeightStructure([[-1], [0], [1]], [[1, -im, 0], [0, 0, 1], [1, im, 0]])
    hw_spaces = [WeightSpace([1], [1, im, 0])]
    return LieAlgebra{Complex{Rational{Int}}, Int}("𝖘𝖔(3)", ch_basis, ws, hw_spaces)
end

# TODO
function LieAlgebra(name::String, size::Int)
    if name == "so" && size == 3
        return so3()
    else
        error("Not implemented")
    end
end

name(alg::LieAlgebra) = alg.name
dim(alg::LieAlgebra) = length(alg.basis.std_basis)
rank(alg::LieAlgebra) = length(alg.basis.cartan)
Base.size(alg::LieAlgebra) = size(alg.basis.std_basis[1], 1)

function Base.show(io::IO, alg::LieAlgebra{F, W}; offset::Int=0) where {F, W}
    println(io, " "^offset, "LieAlgebra $(name(alg))")
    println(io, " "^offset, " number type (or field): $(F)")
    println(io, " "^offset, " weight type: Vector{$(W)}")
    println(io, " "^offset, " dimension: $(dim(alg))")
    print(io, " "^offset, " rank (dimension of Cartan subalgebra): $(rank(alg))")
end

function basis(alg::LieAlgebra; as_matrices::Bool=false)
    if as_matrices
        return alg.basis.std_basis
    end
    coeffs = eye(ComplexF64, dim(alg))
    return [LieAlgebraElem(alg, c) for c in eachcol(coeffs)]
end

cartan_subalgebra(alg::LieAlgebra) = [LieAlgebraElem(alg, coeffs) for coeffs in alg.basis.cartan]
positive_root_elements(
    alg::LieAlgebra
) = [RootElem(alg, coeffs, root) for (coeffs, root) in zip(alg.basis.positive, alg.basis.positive_roots)]
negative_root_elements(
    alg::LieAlgebra
) = [RootElem(alg, coeffs, root) for (coeffs, root) in zip(alg.basis.negative, alg.basis.negative_roots)]
weight_structure(alg::LieAlgebra) = alg.weight_structure
weights(alg::LieAlgebra) = weights(alg.weight_structure)
nweights(alg::LieAlgebra) = nweights(alg.weight_structure)
hw_spaces(alg::LieAlgebra) = alg.hw_spaces


struct SumLieAlgebra <: AbstractReductiveLieAlgebra
    name::String
    algs::Vector{AbstractReductiveLieAlgebra}
end

SumLieAlgebra(
    algs::Vector{AbstractReductiveLieAlgebra}
) = SumLieAlgebra(join([name(alg) for alg in algs], " ⊕ "), algs)

name(alg::SumLieAlgebra) = alg.name
algebras(alg::SumLieAlgebra) = alg.algs
dim(alg::SumLieAlgebra) = sum([dim(a) for a in algebras(alg)])
rank(alg::SumLieAlgebra) = sum([rank(a) for a in algebras(alg)])

function Base.show(io::IO, alg::SumLieAlgebra; offset::Int=0)
    println(io, " "^offset, "SumLieAlgebra $(name(alg))")
    println(io, " "^offset, " dimension: $(dim(alg))")
    print(io, " "^offset, " rank (dimension of Cartan subalgebra): $(rank(alg))")
end

⊕(
    alg₁::AbstractReductiveLieAlgebra,
    alg₂::AbstractReductiveLieAlgebra
) = SumLieAlgebra("$(name(alg₁)) ⊕ $(name(alg₂))", [alg₁, alg₂])

⊕(
    alg₁::SumLieAlgebra,
    alg₂::AbstractReductiveLieAlgebra
) = SumLieAlgebra("$(name(alg₁)) ⊕ $(name(alg₂))", [algebras(alg₁)..., alg₂])

function get_elements(alg::SumLieAlgebra, sym::Symbol)
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

cartan_subalgebra(alg::SumLieAlgebra) = get_elements(alg, :cartan_subalgebra)
positive_root_elements(alg::SumLieAlgebra) = get_elements(alg, :positive_root_elements)
negative_root_elements(alg::SumLieAlgebra) = get_elements(alg, :negative_root_elements)
