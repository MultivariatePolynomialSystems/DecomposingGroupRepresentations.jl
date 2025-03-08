export AbstractSLP, Poly, SumPoly, ProdPoly
export zero_combinations, evaluate, evaluation_matrix
export differentiate

abstract type AbstractSLP{V,M,T} end

Base.zero(::AbstractSLP{V,M,T}) where {V,M,T} = Poly(zero(Polynomial{V,M,T}))
Base.one(::AbstractSLP{V,M,T}) where {V,M,T} = Poly(one(Polynomial{V,M,T}))

struct Poly{V,M,T} <: AbstractSLP{V,M,T}
    poly::Polynomial{V,M,T}
end

Poly(var::Variable) = Poly(polynomial(var))
Base.show(io::IO, p::Poly) = print(io, repr(p.poly))
Base.iszero(p::Poly) = iszero(p.poly)
DynamicPolynomials.isconstant(p::Poly) = isconstant(p.poly)
DynamicPolynomials.variables(p::Poly) = variables(p.poly)
DynamicPolynomials.nvariables(p::Poly) = nvariables(p.poly)
depth(p::Poly) = 1

struct SumPoly{V,M,T} <: AbstractSLP{V,M,T}
    polys::Vector{AbstractSLP{V,M,T}}
end

SumPoly(p::Vector{<:AbstractSLP{V,M,T}}) where {V,M,T} = SumPoly{V,M,T}(p)

Base.show(io::IO, p::SumPoly) = print(io, join([repr(poly) for poly in p.polys], " + "))
DynamicPolynomials.variables(p::SumPoly) = ∪([variables(poly) for poly in p.polys]...)
DynamicPolynomials.nvariables(p::SumPoly) = length(variables(p))
depth(p::SumPoly) = 1 + maximum([depth(poly) for poly in p.polys])

struct ProdPoly{V,M,T,P<:AbstractSLP{V,M,T},S<:AbstractSLP{V,M,T}} <: AbstractSLP{V,M,T}
    p₁::P
    p₂::S
end

Base.show(io::IO, p::ProdPoly) = print(io, "($(repr(p.p₁))) * ($(repr(p.p₂)))")
DynamicPolynomials.variables(p::ProdPoly) = ∪(variables(p.p₁), variables(p.p₂))
DynamicPolynomials.nvariables(p::ProdPoly) = length(variables(p))
depth(p::ProdPoly) = 1 + maximum([depth(p.p₁), depth(p.p₂)])

Base.:+(f::AbstractSLP{V,M,T}, g::AbstractSLP{V,M,T}) where {V,M,T} = iszero(f) ? g : iszero(g) ? f : SumPoly([f, g])
Base.:+(f::Poly, g::Poly) = iszero(f) ? g : iszero(g) ? f : Poly(f.poly + g.poly)
Base.:+(f::SumPoly, g::AbstractSLP) = iszero(g) ? f : SumPoly(vcat(f.polys, [g]))
Base.:+(f::AbstractSLP, g::SumPoly) = iszero(f) ? g : SumPoly(vcat([f], g.polys))
Base.:+(f::SumPoly, g::SumPoly) = SumPoly(vcat(f.polys, g.polys))

function Base.:*(f::AbstractSLP, g::AbstractSLP)
    if iszero(f) || iszero(g)
        return zero(f)
    end
    return ProdPoly(f, g)
end
Base.:*(f::Poly, g::Poly) = iszero(f) || iszero(g) ? zero(f) : isconstant(f) || isconstant(g) ? Poly(f.poly * g.poly) : ProdPoly(f, g)
Base.:*(f::Poly, n::Number) = Poly(f.poly * n)
Base.:*(n::Number, f::Poly) = Poly(f.poly * n)
Base.:*(f::SumPoly, n::Number) = iszero(n) ? zero(f) : SumPoly([poly * n for poly in f.polys])
Base.:*(n::Number, f::SumPoly) = iszero(n) ? zero(f) :  SumPoly([n * poly for poly in f.polys])
Base.:*(f::ProdPoly, n::Number) = iszero(n) ? zero(f) : ProdPoly(f.p₁ * n, f.p₂)
Base.:*(n::Number, f::ProdPoly) = iszero(n) ? zero(f) : ProdPoly(f.p₁ * n, f.p₂)

DynamicPolynomials.differentiate(p::Poly, x::Variable) = Poly(differentiate(p.poly, x))
DynamicPolynomials.differentiate(p::SumPoly, x::Variable) = sum([differentiate(poly, x) for poly in p.polys])
DynamicPolynomials.differentiate(p::ProdPoly, x::Variable) = p.p₁ * differentiate(p.p₂, x) + differentiate(p.p₁, x) * p.p₂

DynamicPolynomials.subs(p::Poly, vars::Vector{<:Variable}, sbs::Vector) = Poly(subs(p.poly, vars => sbs))

evaluate(
    p::AbstractPolynomialLike,
    vars::Vector{<:Variable},
    vals::Vector{T}
) where T = p(vars => vals)
evaluate(
    p::Poly,
    vars::Vector{<:Variable},
    vals::Vector{T}
) where T = evaluate(p.poly, vars, vals)
evaluate(
    p::SumPoly,
    vars::Vector{<:Variable},
    vals::Vector{T}
) where T = sum([evaluate(poly, vars, vals) for poly in p.polys])
evaluate(
    p::ProdPoly,
    vars::Vector{<:Variable},
    vals::Vector{T}
) where T = evaluate(p.p₁, vars, vals) * evaluate(p.p₂, vars, vals)

function evaluation_matrix(exprs::Vector{<:AbstractSLP}, F::DataType; nsamples::Int=1)
    vars = ∪([variables(expr) for expr in exprs]...)
    nexprs, nvars = length(exprs), length(vars)
    sbs = [rand(F, nvars) for _ in 1:nsamples]
    evals = zeros(F, nexprs, nsamples)
    for i in 1:nsamples
        for j in 1:nexprs
            evals[i,j] = evaluate(exprs[j], vars, sbs[i])
        end
    end
    return evals
end

function zero_combinations(exprs::Vector{<:AbstractSLP}, F::DataType; tol::Real=1e-5)
    evals = evaluation_matrix(exprs, F; nsamples=length(exprs))
    N = Matrix(transpose(nullspace(evals; atol=tol)))
    rref!(N)
    sparsify!(N, tol)
    return eachrow(N)
end

Base.isapprox(f::Poly, g::Poly; kwargs...) = isapprox(f.poly, g.poly; kwargs...)

function Base.isapprox(f::AbstractSLP, g::AbstractSLP; nsamples::Int=1, kwargs...)
    vars = variables(f) ∪ variables(g)
    # isempty(vars) && return isapprox(convert(ComplexF64, f), convert(ComplexF64, g); kwargs...)
    sbs = [rand(ComplexF64, length(vars)) for _ in 1:nsamples]
    evals_f = [evaluate(f, vars, sb) for sb in sbs]
    evals_g = [evaluate(g, vars, sb) for sb in sbs]
    return isapprox(evals_f, evals_g; kwargs...)
end