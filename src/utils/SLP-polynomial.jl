export MyAbstractPoly, MyPoly, SumPoly, ProdPoly

abstract type MyAbstractPoly end

struct MyPoly{T <: Polynomial} <: MyAbstractPoly
    poly::T
end

function Base.show(io::IO, p::MyPoly)
    println(io, "MyPoly")
    print(io, repr(p.poly))
end

struct SumPoly <: MyAbstractPoly
    polys::Vector{MyAbstractPoly}
end

function Base.show(io::IO, p::SumPoly)
    println(io, "SumPoly")
    for poly in p.polys
        print(io, "poly: ")
        show(io, poly)
        println(io)
    end
end

struct ProdPoly{T <: MyAbstractPoly, S <: MyAbstractPoly} <: MyAbstractPoly
    p₁::T
    p₂::S
end

function Base.show(io::IO, p::ProdPoly)
    println(io, "ProdPoly")
    print(io, "p₁: ")
    show(io, p.p₁)
    println(io)
    print(io, "p₂: ")
    show(io, p.p₂)
end

Base.:+(f::MyAbstractPoly, g::MyAbstractPoly) = SumPoly([f, g])
Base.:+(f::SumPoly, g::MyAbstractPoly) = SumPoly(vcat(f.polys, [g]))
Base.:+(f::MyAbstractPoly, g::SumPoly) = SumPoly(vcat([f], g.polys))
Base.:+(f::SumPoly, g::SumPoly) = SumPoly(vcat(f.polys, g.polys))

Base.:*(f::MyAbstractPoly, g::MyAbstractPoly) = ProdPoly(f, g)

DynamicPolynomials.differentiate(p::MyPoly, x::Variable) = MyPoly(differentiate(p.poly, x))
DynamicPolynomials.differentiate(p::SumPoly, x::Variable) = SumPoly([differentiate(poly, x) for poly in p.polys])
DynamicPolynomials.differentiate(p::ProdPoly, x::Variable) = p.p₁ * differentiate(p.p₂, x) + differentiate(p.p₁, x) * p.p₂

evaluate(
    p::MyPoly,
    vars::NTuple{N, <:Variable},
    vals::NTuple{N, <:Number}
) where N = p.poly(vars => vals)
evaluate(
    p::SumPoly,
    vars::NTuple{N, <:Variable},
    vals::NTuple{N, <:Number}
) where N = sum([evaluate(f, vars, vals) for f in p.polys])
evaluate(
    p::ProdPoly,
    vars::NTuple{N, <:Variable},
    vals::NTuple{N, <:Number}
) where N = evaluate(p.p₁, vars, vals) * evaluate(p.p₂, vars, vals)
