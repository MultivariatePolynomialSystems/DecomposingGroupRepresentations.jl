

a2p(M::AbstractMatrix{<:Number}) = [M; ones(eltype(M), 1, size(M, 2))]
p2a(M::AbstractMatrix{<:Number}) = (M./M[end:end,:])[1:end-1,:]

M2VV(M::AbstractMatrix; copy::Bool=true) = copy ? [M[:,i] for i in axes(M, 2)] : [view(M,:,i) for i in axes(M,2)]

V2M(v::AbstractVector) = reshape(v, length(v), 1) # TODO: do I need this? I can access a Vector with [i,1] too...
M2V(M::AbstractMatrix) = M[:,1] # TODO: same here, I can access elements of 1-column matrix M by M[i] too...
