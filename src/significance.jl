#####################################################
# Fourier surrogates
#####################################################

# Shift each frequency content by a random phase.
function generate_fourier_surrogate(x::Vector{T}) where {T<:Real}
    F = rfft(x)
    return irfft(F .* exp.(2 * π * im .* rand(length(F))), length(x))
end

"""
    generate_stacked_fourier_surrogates(X::AbstractArray, ns::Int)

Generate ns Fourier surrogates for each line of X.
If X is a CuArray, computation takes place in GPU.
Aformentionned speedup particularly significant if X big or ns large.
"""
function generate_stacked_fourier_surrogates(x::Vector{T}, ns::Int) where {T<:Real}
    S = zeros(T, ns, length(x))
    for i in axes(S, 1)
        S[i, :] = generate_fourier_surrogate(x)
    end
    return S
end

function generate_stacked_fourier_surrogates(X::Matrix{T}, ns::Int) where {T<:Real}
    nx, nt = size(X)
    F = repeat(rfft(X, 2), inner = (ns, 1))
    stacked_surrogates = irfft(
        F .* exp.(T(2 * π) * im .* rand(T, nx * ns, size(F, 2))),
        nt,
        2,
    )
    return stacked_surrogates
end

function generate_stacked_fourier_surrogates(X::CuArray{T,2}, ns::Int) where {T<:Real}
    nx, nt = size(X)
    F = repeat(CUDA.CUFFT.rfft(X, 2), inner = (ns, 1))
    stacked_surrogates = CUDA.CUFFT.irfft(
        F .* exp.(T(2 * π) * im .* CUDA.rand(T, nx * ns, size(F, 2))),
        nt,
        2,
    )
    return stacked_surrogates
end

#####################################################
# Percentile significance 
#####################################################
"""
    percentile_significance(ref_stat::AbstractArray, sur_stat::Matrix{T}, ns::Int, nx::Int)

Get percentile significance of reference (original time series) versus surrogate statistics.
"""
function percentile_significance(
    ref_stat::Matrix{T},
    sur_stat::Matrix{T},
    ns::Int,
    nx::Int,
) where {T<:Real}

    rep_ref = repeat(ref_stat, inner = (ns, 1))
    return (kron(I(nx), ones(T, ns)') * (sur_stat .- rep_ref .< T(0))) ./ T(ns)
end

function percentile_significance(
    ref_stat::CuArray{T,2},
    sur_stat::CuArray{T,2},
    ns::Int,
    nx::Int,
) where {T<:Real}

    rep_ref = repeat(ref_stat, inner = (ns, 1))
    return (CuArray(kron(I(nx), ones(T, ns)')) * (sur_stat .- rep_ref .< T(0))) ./ T(ns)
end

#####################################################
# Predict transition
#####################################################

# TODO insert tolerance wrt lag
function count_positive_indicators(
    indicator_trend_significance3D::Union{Array{T, 3}, CuArray{T, 3}};
    plevel=T(0.95),
    nindicators::Int=size(indicator_trend_significance3D, 3),
    threshold::Bool=false,
) where {T}
    Σ = sum_significant_indicators(indicator_trend_significance3D, plevel)
    prediction = Σ ./ nindicators
    if threshold
        prediction = threshold_indicator_significance(prediction)
    end
    return prediction
end

function stack_indicators(indicator_list::Any)
    return cat(indicator_list..., dims=3)
end

function sum_percentiles(P::Union{Array{T, 3}, CuArray{T,3}}) where {T}
    return reduce( +, P, dims=3 )[:,:,1]
end

function sum_significant_indicators(P::Union{Array{T, 3}, CuArray{T,3}}, plevel::T) where {T}
    return reduce( +, T.(P .> plevel), dims=3 )[:,:,1]
end

function threshold_indicator_significance(S::Union{Matrix{T}, CuArray{T,2}}) where {T}
    return isapprox.(S, 1)
end