export rand_rotation

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
end

function sparsify!(M::AbstractMatrix{<:Number}, tol::Real; digits::Integer=0)
    for r in eachrow(M)
        sparsify!(r, tol; digits=digits)
    end
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