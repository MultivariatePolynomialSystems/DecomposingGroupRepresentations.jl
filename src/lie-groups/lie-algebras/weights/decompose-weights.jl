export sym_weight_dict

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

# TODO: works only for so(3)
hw2all(hw::Weight, ::LieAlgebra) = [Weight([j]) for j in -hw[1]:hw[1]]
hw2all(hw::Weight, ::ScalingLieAlgebra) = [hw]

function hw2all(hw::Weight, A::SumLieAlgebra)
    ws = split_weight(hw, A)
    all_ws = [hw2all(w, alg) for (w, alg) in zip(ws, algebras(A))]
    all_ws = product(all_ws...)
    return [vcat(collect(ws)...) for ws in all_ws][:]
end

function sym_weight_dict(hw::T, A::AbstractLieAlgebra, p::Int) where T<:Weight
    ws = hw2all(hw, A)
    combs = multiexponents(; degree=p, nvars=length(ws))
    ws_dict = Dict{T, Int}()
    for comb in combs
        w = sum([val*ws[ind] for (ind, val) in zip(comb.nzind, comb.nzval)])
        ws_dict[w] = get(ws_dict, w, 0) + 1
    end
    return p*hw, ws_dict
end

# # TODO: works only for so(3)
# child_weights(hw::T, ::LieAlgebra) where T<:Weight = iszero(hw[1]) ? T[] : [Weight([hw[1]-1])]
# child_weights(::T, ::ScalingLieAlgebra) where T<:Weight = T[]
# function child_weights(hw::Weight, A::SumLieAlgebra)
#     ws = split_weight(hw, A)
#     all_child_ws = [child_weights(w, alg) for (w, alg) in zip(ws, algebras(A))]
#     full_child_ws = [[copy(ws) for _ in child_ws] for child_ws in all_child_ws]
#     for (i, child_ws) in enumerate(all_child_ws)
#         for (j, w) in enumerate(child_ws)
#             full_child_ws[i][j][i] = w
#         end
#     end
#     return [vcat(ch_ws...) for ch_ws in flatten(full_child_ws)]
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
    decomp_hws::Union{Set{W}, Vector{W}}
) where W<:Weight
    if ws_dict[hw] != 0
        orbit_ws = hw2all(hw, A)
        while ws_dict[hw] > 0
            push!(decomp_hws, hw)
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

function sym(hw::T, A::AbstractLieAlgebra, p::Int; with_muls::Bool=false) where T<:Weight
    sym_hw, ws_dict = sym_weight_dict(hw, A, p)
    decomp_hws = with_muls ? T[] : Set{T}()
    sym!(sym_hw, A, ws_dict, decomp_hws)
    return collect(decomp_hws)
end
