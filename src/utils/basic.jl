export rand_rotation, sparsify!
export rref, monomials, coeffs_matrix, polynomials
export num_mons
export superscript, subscript
export multiexponents, div_by_smallest_coeff
export eye, a2p, p2a, rand_unit

a2p(M::AbstractMatrix{<:Number}) = [M; ones(eltype(M), 1, size(M, 2))]
p2a(M::AbstractMatrix{<:Number}) = (M./M[end:end,:])[1:end-1,:]

M2VV(M::AbstractMatrix; copy::Bool=true) = copy ? [M[:,i] for i in axes(M, 2)] : [view(M,:,i) for i in axes(M,2)]

V2M(v::AbstractVector) = reshape(v, length(v), 1) # TODO: do I need this? I can access a Vector with [i,1] too...
M2V(M::AbstractMatrix) = M[:,1] # TODO: same here, I can access elements of 1-column matrix M by M[i] too...

Base.rand(::Type{Rational{T}}) where T <: Integer = rand(T(-10):T(10)) // rand(T(-10):T(10))
Base.rand(::Type{Complex{Rational{T}}}) where T <: Integer = rand(Rational{T}) + rand(Rational{T})*im
Base.rand(T::DataType, n::Int) = [rand(T) for _ in 1:n]

eye(n::Integer) = eye(Float64, n)
eye(T::Type, n::Integer) = Matrix{T}(I(n))

function q2sR(q)
    a,b,c,d = q
    return [a^2+b^2-c^2-d^2 2*(b*c-a*d) 2*(b*d+a*c);
            2*(b*c+a*d) a^2-b^2+c^2-d^2 2*(c*d-a*b);
            2*(b*d-a*c) 2*(c*d+a*b) a^2-b^2-c^2+d^2]
end

function q2R(q)
    a,b,c,d = q
    return 1/(a^2+b^2+c^2+d^2)*q2sR(q)
end

rand_rotation(T::DataType) = q2R(rand(T, 4))

function num_mons(n::Integer, d::Integer; upto::Bool=false)
    upto && return n > 0 ? binomial(Int(n + d), Int(d)) : 0
    return n > 0 ? binomial(Int(n - 1 + d), Int(d)) : 0
end

# TODO: test this
function sparsify!(v::AbstractVector{<:Number}, tol::Real; digits::Integer=0)
    for j in eachindex(v)
        if abs(imag(v[j])) < tol
            v[j] = real(v[j])
        elseif abs(round(imag(v[j]); digits=digits) - imag(v[j])) < tol
            v[j] = real(v[j]) + round(imag(v[j]); digits=digits)*im
        end
        if abs(real(v[j])) < tol
            v[j] = imag(v[j])*im
        elseif abs(round(real(v[j]); digits=digits) - real(v[j])) < tol
            v[j] = round(real(v[j]); digits=digits) + imag(v[j])*im
        end
    end
    return v
end

function sparsify!(M::AbstractMatrix{<:Number}, tol::Real; digits::Integer=0)
    for r in eachrow(M)
        sparsify!(r, tol; digits=digits)
    end
    return M
end

function simplify_numbers(v::AbstractVector{<:Number})
    v = Vector{Number}(v)
    for (i, vᵢ) in enumerate(v)
        try
            v[i] = Integer(vᵢ)
        catch
            try
                v[i] = Real(vᵢ)
            catch
                try
                    v[i] = Complex{Integer}(vᵢ)
                catch
                end
            end
        end
    end
    return v
end

function div_by_lowest_magnitude(v::AbstractVector{<:Number}, tol::Float64)
    j = 1
    while norm(v[j]) < tol
        j += 1
        if j > length(v)
            return v
        end
    end
    a = v[j]
    for vᵢ in v
        if tol < norm(vᵢ) < norm(a)
            a = vᵢ
        end
    end
    new_v = v./a
    sparsify!(new_v, tol)
    return new_v
end

function subscript(n::Integer)::String
    c = n < 0 ? [Char(0x208B)] : []
    for d in reverse(digits(abs(n)))
        push!(c, Char(0x2080+d))
    end
    return join(c)
end

function superscript(n::Integer)::String
    c = n < 0 ? [Char(0x207B)] : []
    for d in reverse(digits(abs(n)))
        if d == 0 push!(c, Char(0x2070)) end
        if d == 1 push!(c, Char(0x00B9)) end
        if d == 2 push!(c, Char(0x00B2)) end
        if d == 3 push!(c, Char(0x00B3)) end
        if d > 3 push!(c, Char(0x2070+d)) end
    end
    return join(c)
end

DynamicPolynomials.monomials(
    F::Vector{<:AbstractPolynomialLike};
    as_monvec::Bool=true
) = as_monvec ? MonomialVector(∪([monomials(f) for f in F]...)) : ∪([monomials(f) for f in F]...)
DynamicPolynomials.monomials(Fs::Vector{<:AbstractPolynomialLike}...) = MonomialVector(∪([monomials(F) for F in Fs]...))

# column corresponds to a polynomial
function coeffs_matrix(
    F::Vector{<:AbstractPolynomialLike{T}},
    mons::AbstractVector{M}
) where {T, M<:Monomial}
    C = zeros(T, length(mons), length(F))
    for (i, f) in enumerate(F)
        d = Dict{M,T}(zip(monomials(f), coefficients(f)))
        for (j, mon) in enumerate(mons)
            C[j, i] = get(d, mon, zero(T))
        end
    end
    return C
end

# polynomials(
#     M::AbstractMatrix,
#     mons::Vector{<:Monomial}
# ) = [sum(m .* mons) for m in eachcol(M)]

# polynomials(
#     C::AbstractVector{<:AbstractVector},
#     mons::Vector{<:Monomial}
# ) = [sum(c .* mons) for c in C]

polynomials(
    C::AbstractVector{<:AbstractVector{T}},
    mons::MonomialVector{V,M}
) where {V,M,T} = [Polynomial{V,M,T}(c, copy(mons)) for c in C]

function rref(
    F::Vector{<:AbstractPolynomial{T}};
    tol::Real=1e-5
) where T
    mons = monomials(F)
    M = Matrix(transpose(coeffs_matrix(F, mons)))
    rref!(M, tol)
    sparsify!(M, tol)
    N = filter(row -> !iszero(row), eachrow(M))
    return polynomials(N, mons)
end

function zero_combinations(
    F::Vector{<:AbstractPolynomial{T}};
    tol::Real=1e-5,
    logging::Bool=true,
    square_up::Bool=true
) where T
    mons = monomials(F; as_monvec=false)
    M = coeffs_matrix(F, mons)
    # @assert size(M, 1) ≥ size(M, 2)
    if size(M, 1) > size(M, 2) && square_up
        logging && print(crayon"#f4d03f", "Squaring up the $(size(M)) coefficients matrix...\n")
        M = rand(T, size(M, 2), size(M, 1)) * M
    end
    logging && print(crayon"#f4d03f", "Computing nullspace of $(size(M)) matrix...\n")
    N = Matrix(transpose(nullspace(M; atol=tol)))
    logging && print(crayon"#f4d03f", "Nullspace is $(size(N, 1))-dimensional\n")
    sparsify!(N, tol)
    rref!(N, tol)
    display(N)
    for n in N
        if norm(n) > 1e3
            display(N)
            error("Nullspace contains large coefficients")
        end
    end
    sparsify!(N, tol)
    return eachrow(N)
end

# Generates a list of multiexponents of degree @degree in @nvars variables
function multiexponents(degree::Tv, nvars::Ti) where {Tv<:Integer,Ti<:Integer}
    mexps = [spzeros(Tv, Ti, nvars) for _ in 1:num_mons(nvars, degree)]
    iszero(degree) && return mexps
    k = 1
    for n in 1:nvars
        for part::Vector{Tv} in partitions(degree, n)
            for vals in multiset_permutations(part, n)
                for inds in combinations(Ti.(1:nvars), n)
                    mexps[k][inds] = vals
                    k += 1
                end
            end
        end
    end
    return mexps
end

function multiexponents(; degree::Tv, nvars::Ti, upto::Bool=false) where {Tv<:Integer,Ti<:Integer}
    !upto && return multiexponents(degree, nvars)
    mexps = [spzeros(Tv, Ti, nvars)]
    for d::Tv in 1:degree
        append!(mexps, multiexponents(d, nvars))
    end
    return mexps
end

function div_by_smallest_coeff(f::Polynomial; tol::Real=1e-5)
    cs = coefficients(f)
    sparsify!(cs, tol)
    nzcs = filter(x -> abs(x) > tol, cs)
    c = minimum(abs, nzcs)
    cs = cs/c
    sparsify!(cs, tol)
    return Polynomial(cs, monomials(f))
end

function vec_subscript(v::AbstractVector{<:Integer})
    return join([subscript(vᵢ) for vᵢ in v], "ˏ")
end

function rand_unit(dims...)
    r = rand(ComplexF64, dims...)
    for i in eachindex(r)
        r[i] = r[i]/abs(r[i])
    end
    return r
end