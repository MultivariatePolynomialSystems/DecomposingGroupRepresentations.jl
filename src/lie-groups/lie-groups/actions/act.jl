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

act(
    X::LieAlgebraElem{F, T},
    a::MatrixGroupAction{Lie, F},
    f::AbstractPolynomialLike
) where {F, T<:AbstractLieAlgebra{F}} = sum([sum(differentiate(f, vars) .* -matrix(X)*vars) for vars in action_vectors(a)])

act(
    X::LieAlgebraElem{F, ScalingLieAlgebra{F}},
    a::ScalingLieGroupAction{F},
    f::AbstractPolynomialLike
) where F = act(X, MatrixGroupAction(a), f)

act(
    X::SumLieAlgebraElem{F},
    a::DirectProductGroupAction{Lie, F},
    f::AbstractPolynomialLike
) where F = sum([act(Xᵢ, aᵢ, f) for (Xᵢ, aᵢ) in zip(elements(X), lie_actions(a))])

act(
    elem::RootElem,
    a::AbstractGroupAction,
    wv::WeightVector
) = WeightVector(root(elem)+weight(wv), act(elem.elem, a, vector(wv)))