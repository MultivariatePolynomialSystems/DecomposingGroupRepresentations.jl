export ScalingLieAlgebra,
    LieAlgebra,
    weight_structure,
    cartan_subalgebra,
    positive_root_elements,
    negative_root_elements,
    SumLieAlgebra,
    algebras


struct ScalingLieAlgebra{F} <: AbstractLieAlgebra{F}
    name::String
    exps::SparseMatrixCSC{Int, Int} # every row is a vector u = [u₁,...,uₖ] which acts on vars by λᵘ
end

ScalingLieAlgebra{F}(size::Int) where F = ScalingLieAlgebra{F}("ℂ", sparse(ones(Int, 1, size)))

name(alg::ScalingLieAlgebra) = alg.name
dim(alg::ScalingLieAlgebra) = size(alg.exps, 1)
exponents(alg::ScalingLieAlgebra) = alg.exps
rank(alg::ScalingLieAlgebra) = dim(alg)
Base.size(alg::ScalingLieAlgebra) = size(alg.exps, 2)
weight_type(::ScalingLieAlgebra) = Weight{Int}

function show_basis(io::IO, alg::ScalingLieAlgebra; offset::Int=0)
    for i in 1:dim(alg)
        print(io, " "^offset, "[", join(alg.exps[i, :], ", "), "]")
        i < dim(alg) && print(io, '\n')
    end
end

function Base.show(io::IO, alg::ScalingLieAlgebra{F}; offset::Int=0) where F
    println(io, " "^offset, "ScalingLieAlgebra $(name(alg))")
    println(io, " "^offset, " number type (or field): $(F)")
    println(io, " "^offset, " dimension: $(dim(alg))")
    println(io, " "^offset, " basis (diagonal matrices):")
    show_basis(io, alg, offset=offset+2)
end

function basis(alg::ScalingLieAlgebra{F}; as_matrices::Bool=false) where F
    as_matrices && return [Diagonal(F.(e)) for e in eachrow(alg.exps)]
    coeffs = eye(F, dim(alg))
    return [LieAlgebraElem(alg, c) for c in eachcol(coeffs)]
end

cartan_subalgebra(alg::ScalingLieAlgebra) = basis(alg)
positive_root_elements(::ScalingLieAlgebra) = LieAlgebraElem{ScalingLieAlgebra}[]
negative_root_elements(::ScalingLieAlgebra) = LieAlgebraElem{ScalingLieAlgebra}[]

struct Root
    root::Vector{Int}
end

Base.:+(r::Root, w::Weight) = Weight(r.root + w.weight)
Base.:+(w::Weight, r::Root) = Weight(w.weight + r.root)
Base.:+(r₁::Root, r₂::Root) = Root(r₁.root + r₂.root)
Base.convert(::Type{Root}, v::Vector{Int}) = Root(v)
Base.convert(::Type{Vector{Int}}, r::Root) = r.root
to_vector(r::Root) = r.root

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


struct LieAlgebra{F, W<:Weight} <: AbstractLieAlgebra{F}
    name::String
    basis::ChevalleyBasis{F}
    weight_structure::WeightStructure{MatrixVectorSpace{F}, W}
    hw_spaces::Vector{WeightSpace{MatrixVectorSpace{F}, W}} # TODO: change to WeightStructure?
end

function so3(field_type::DataType, weight_type::DataType)
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
    return LieAlgebra{field_type, Weight{weight_type}}("𝖘𝖔(3)", ch_basis, ws, hw_spaces)
end
so3() = so3(ComplexF64, Int)

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
weight_type(::LieAlgebra{F, W}) where {F, W} = W

function Base.show(io::IO, alg::LieAlgebra{F, Weight{W}}; offset::Int=0) where {F, W}
    println(io, " "^offset, "LieAlgebra $(name(alg))")
    println(io, " "^offset, " number type (or field): $(F)")
    println(io, " "^offset, " weight type: $(W)")
    println(io, " "^offset, " dimension: $(dim(alg))")
    print(io, " "^offset, " rank (dimension of Cartan subalgebra): $(rank(alg))")
end

function basis(alg::LieAlgebra{F}; as_matrices::Bool=false) where F
    if as_matrices
        return alg.basis.std_basis
    end
    coeffs = eye(F, dim(alg))
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


struct SumLieAlgebra{F, T<:AbstractLieAlgebra{F}} <: AbstractLieAlgebra{F}
    name::String
    algs::Vector{T}
end

SumLieAlgebra(
    algs::Vector{<:AbstractLieAlgebra{F}}
) where F = SumLieAlgebra(join([name(alg) for alg in algs], " ⊕ "), algs)

name(alg::SumLieAlgebra) = alg.name
algebras(alg::SumLieAlgebra) = alg.algs
dim(alg::SumLieAlgebra) = sum([dim(a) for a in algebras(alg)])
rank(alg::SumLieAlgebra) = sum([rank(a) for a in algebras(alg)])
weight_type(alg::SumLieAlgebra) = promote_type([weight_type(a) for a in algebras(alg)]...)

function Base.show(io::IO, alg::SumLieAlgebra{F}; offset::Int=0) where F
    println(io, " "^offset, "SumLieAlgebra $(name(alg))")
    println(io, " "^offset, " number type (or field): $(F)")
    println(io, " "^offset, " dimension: $(dim(alg))")
    print(io, " "^offset, " rank (dimension of Cartan subalgebra): $(rank(alg))")
end

⊕(
    alg₁::AbstractLieAlgebra{F},
    alg₂::AbstractLieAlgebra{F}
) where F = SumLieAlgebra("$(name(alg₁)) ⊕ $(name(alg₂))", [alg₁, alg₂])

⊕(
    alg₁::SumLieAlgebra{F},
    alg₂::AbstractLieAlgebra{F}
) where F = SumLieAlgebra("$(name(alg₁)) ⊕ $(name(alg₂))", [algebras(alg₁)..., alg₂])

function get_elements(alg::SumLieAlgebra, sym::Symbol)
    if sym == :positive_root_elements || sym == :negative_root_elements
        elems = RootElem[]
    else
        elems = SumLieAlgebraElem[]
    end
    for (i, a) in enumerate(algebras(alg))
        a_elems = eval(sym)(a)
        alg_elems = [zero(alg) for _ in a_elems]
        if sym == :positive_root_elements || sym == :negative_root_elements
            roots = [[zeros(Int, rank(a)) for a in algebras(alg)] for _ in a_elems]
        end
        for (j, elem) in enumerate(a_elems)
            if sym == :positive_root_elements || sym == :negative_root_elements
                alg_elems[j][i] = element(elem)
                roots[j][i] = root(elem)
            else
                alg_elems[j][i] = elem
            end
        end
        if sym == :positive_root_elements || sym == :negative_root_elements
            append!(elems, [RootElem(elem, vcat(root...)) for (elem, root) in zip(alg_elems, roots)])
        else
            append!(elems, alg_elems)
        end
    end
    return elems
end

cartan_subalgebra(alg::SumLieAlgebra) = get_elements(alg, :cartan_subalgebra)
positive_root_elements(alg::SumLieAlgebra) = get_elements(alg, :positive_root_elements)
negative_root_elements(alg::SumLieAlgebra) = get_elements(alg, :negative_root_elements)

zero_weight(alg::AbstractLieAlgebra) = zero(weight_type(alg), rank(alg))
positive_roots(alg::AbstractLieAlgebra) = [root(pre) for pre in positive_root_elements(alg)]
negative_roots(alg::AbstractLieAlgebra) = [root(nre) for nre in negative_root_elements(alg)]