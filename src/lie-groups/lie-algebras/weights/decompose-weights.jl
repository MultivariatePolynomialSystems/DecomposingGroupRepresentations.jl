export sym_weight_muls_dict

# TODO: works only for so(3)
function ⊗(hw₁::Weight, hw₂::Weight, ::LieAlgebra)
    @assert length(hw₁) == length(hw₂) == 1
    m, n = hw₁[1], hw₂[1]
    m, n = m > n ? (m, n) : (n, m)
    return [Weight([m+n-j]) for j in 0:2*n]
end

⊗(hw₁::Weight, hw₂::Weight, ::ScalingLieAlgebra) = [hw₁ + hw₂]


function split_weight(w::T, A::SumLieAlgebra) where T<:Weight
    inds = cumsum([rank(alg) for alg in algebras(A)])
    st_ind = 1
    ws = T[]
    for end_ind in inds
        push!(ws, Weight(w[st_ind:end_ind]))
        st_ind = end_ind + 1
    end
    return ws
end

function ⊗(hw₁::Weight, hw₂::Weight, A::SumLieAlgebra)
    ws₁ = split_weight(hw₁, A)
    ws₂ = split_weight(hw₂, A)
    sep_ws = [⊗(w₁, w₂, alg) for (w₁, w₂, alg) in zip(ws₁, ws₂, algebras(A))]
    all_ws = product(sep_ws...)
    return [vcat(collect(ws)...) for ws in all_ws][:]
end

function tensor(hw₁::Weight, hw₂::Weight, A::AbstractLieAlgebra)
    ws = ⊗(hw₁, hw₂, A)
    return Dict(zip(ws, ones(Int, length(ws)))) # TODO: Each weight has multiplicity 1
end

# TODO: works only for so(3)
hw2all(hw::Weight, ::LieAlgebra) = [Weight([j]) for j in -hw[1]:hw[1]]
hw2all(hw::Weight, ::ScalingLieAlgebra) = [hw]

function hw2all(hw::Weight, A::SumLieAlgebra)
    ws = split_weight(hw, A)
    all_ws = [hw2all(w, alg) for (w, alg) in zip(ws, algebras(A))]
    all_ws = product(all_ws...)
    return [vcat(collect(ws)...) for ws in all_ws][:]
end

# Returns the dictionary of weights of Symᵖ(V(hw)) and the highest weight of Symᵖ(V(hw))
# The values of the dictionary are the multiplicities of the corresponding weights
function sym_weight_muls_dict(hw::T, A::AbstractLieAlgebra, p::Int) where T<:Weight
    ws = hw2all(hw, A)
    combs = multiexponents(; degree=p, nvars=length(ws))
    ws_dict = Dict{T, Int}()
    for comb in combs
        w = sum([val*ws[ind] for (ind, val) in zip(comb.nzind, comb.nzval)])
        ws_dict[w] = get(ws_dict, w, 0) + 1
    end
    return p*hw, ws_dict
end

# function tensor_weight_muls_dict(hw₁::T, hw₂::T, A::AbstractLieAlgebra) where T<:Weight
#     ws₁ = hw2all(hw₁, A)
#     ws₂ = hw2all(hw₂, A)
#     ws_dict = Dict{T, Int}()
#     ws_prod = product(ws₁, ws₂)
#     for ws in ws_prod
#         w = ws[1] + ws[2]
#         ws_dict[w] = get(ws_dict, w, 0) + 1
#     end
#     return ws_dict
# end

function some_highest_weight(
    ws_dict::Dict{W, Int},
    A::AbstractLieAlgebra
) where W<:Weight
    w = first(keys(ws_dict))
    while true
        exists = false
        for α in positive_roots(A)
            if haskey(ws_dict, w + α)
                w = w + α
                exists = true
                break
            end
        end
        !exists && return w
    end
end

function sym!(
    hw::W,
    A::AbstractLieAlgebra,
    ws_dict::Dict{W, Int},
    decomp_hws::Dict{W, Int}
) where W<:Weight
    if ws_dict[hw] != 0
        orbit_ws = hw2all(hw, A)
        while ws_dict[hw] > 0
            decomp_hws[hw] = get(decomp_hws, hw, 0) + 1
            [ws_dict[w] -= 1 for w in orbit_ws]
        end
        for w in collect(keys(ws_dict))
            ws_dict[w] == 0 && delete!(ws_dict, w)
        end
    end
    isempty(ws_dict) && return decomp_hws
    shw = some_highest_weight(ws_dict, A)
    sym!(shw, A, ws_dict, decomp_hws)
end

function sym(hw::W, A::AbstractLieAlgebra, p::Int) where W<:Weight
    sym_hw, ws_dict = sym_weight_muls_dict(hw, A, p)
    decomp_hws = Dict{W, Int}() # for the weight key gives its multiplicity
    sym!(sym_hw, A, ws_dict, decomp_hws)
    return decomp_hws
end
