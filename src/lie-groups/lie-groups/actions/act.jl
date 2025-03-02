function act(
    g::GroupElem{S},
    a::MatrixGroupAction{T, F, S},
    f::AbstractPolynomialLike
) where {T<:GroupType, F, S<:AbstractGroup{T, F}}
    sbs = vcat([matrix(g)*v for v in action_vectors(a)]...)
    vars = vcat(action_vectors(a)...)
    return subs(f, vars => sbs)
end

act(
    g::GroupElem{T},
    a::ScalingLieGroupAction{F, T},
    f::AbstractPolynomialLike
) where {F, T <: ScalingLieGroup{F}} = act(g, MatrixGroupAction(a), f)

function act(
    g::DirectProductGroupElem{S},
    a::DirectProductGroupAction{T, F, S},
    f::AbstractPolynomialLike
) where {T<:GroupType, F, S<:AbstractDirectProductGroup{T, F}}
    fₐ = f
    for (gᵢ, aᵢ) in zip(elements(g), actions(a))
        fₐ = act(gᵢ, aᵢ, fₐ)
    end
    return fₐ
end

function act(
    X::LieAlgebraElem{F, T},
    a::MatrixGroupAction{Lie, F},
    f::AbstractPolynomialLike
) where {F, T<:AbstractLieAlgebra{F}}
    res = zero(f)
    for vars in action_vectors(a)
        mv = matrix(X)*vars
        for (i, p) in enumerate(mv)
            if !iszero(p)
                res += p * differentiate(f, vars[i])
            end
        end
    end
    return res
end

# act(
#     X::LieAlgebraElem{F, T},
#     a::MatrixGroupAction{Lie, F},
#     f::AbstractPolynomialLike
# ) where {F, T<:AbstractLieAlgebra{F}} = sum([sum(differentiate(f, vars) .* -matrix(X)*vars) for vars in action_vectors(a)])

act(
    X::LieAlgebraElem{F, ScalingLieAlgebra{F}},
    a::ScalingLieGroupAction{F},
    f::AbstractPolynomialLike
) where F = act(X, MatrixGroupAction(a), f)

function act(
    X::SumLieAlgebraElem{F},
    a::DirectProductGroupAction{Lie, F},
    f::AbstractPolynomialLike
) where F
    res = zero(f)
    for (Xᵢ, aᵢ) in zip(elements(X), lie_actions(a))
        if !iszero(Xᵢ)
            res += act(Xᵢ, aᵢ, f)
        end
    end
    return res
end

# act(
#     X::SumLieAlgebraElem{F},
#     a::DirectProductGroupAction{Lie, F},
#     f::AbstractPolynomialLike
# ) where F = sum([act(Xᵢ, aᵢ, f) for (Xᵢ, aᵢ) in zip(elements(X), lie_actions(a))])

act(
    X::RootElem,
    a::AbstractGroupAction,
    f::AbstractPolynomialLike
) = act(element(X), a, f)

act(
    X::RootElem,
    a::AbstractGroupAction,
    wv::WeightVector
) = WeightVector(root(X)+weight(wv), act(X, a, vector(wv)))