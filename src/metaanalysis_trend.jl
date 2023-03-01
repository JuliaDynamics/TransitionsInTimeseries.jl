#####################################################
# Trend metrics
#####################################################
"""

    ridge(t, x; lambda = 0)

Perform ridge regression of `x` over `t` with regularization parameter `lambda`.
Return vector containing slope and offset.
If `lambda = 0`, linear regression is recovered (default case).
For more information, visit: https://en.wikipedia.org/wiki/Ridge_regression
"""
function ridge(
    t::AbstractVector{T},
    x::AbstractVector{T};
    lambda::T = T(0),
) where {T<:Real}
    M = precompute_ridge(t, lambda = lambda)
    return M * x
end

function ridge(
    x::AbstractVector{T};
    lambda::T = T(0),
) where {T<:Real}
    return ridge(T.(eachindex(x)), x, lambda=lambda)
end

"""

    ridge_slope(t, x; lambda)

Extract the slope of ridge regression of `x` over `t`.
"""
ridge_slope(args...; kwargs...) = ridge(args...; kwargs...)[1]

"""

    precomputed_ridge_slope(p; lambda)

Return a function that computes the slope obtained by ridge regression of an input `x`.
For performance, some terms are pre-computed based on `p::MetaAnalysisParameters`.
"""
function precomputed_ridge_slope(
    p::MetaAnalysisParameters;
    lambda::Real = 0.0,
)
    m = precompute_ridge_slope(1.0:p.wv_evolution_width+1, lambda = lambda)
    curry(f, y) = x -> f(x, y)
    return curry(precomputed_ridge_slope, m)
end

function precomputed_ridge_slope(
    x::AbstractVector{T},
    m::AbstractVector{T},
) where {T<:AbstractFloat} 
    return m' * x
end

"""

    precompute_ridge(t, lambda = 0)

Precompute the matrix arising in the ridge regression. Particularly suited
if ridge regression is to be computed many times on evenly spaced time.
"""
function precompute_ridge(
    t::AbstractVector{T};
    lambda::T = T(0),
) where {T<:Real}
    TT = hcat(t, ones(T, length(t)))'
    return inv(TT * TT' + lambda .* LinearAlgebra.I(2) ) * TT
end

"""

    precompute_ridge_slope(t, lambda = 0)

Precompute the vector arising in the ridge regression of the slope. Particularly suited
if ridge regression is to be computed many times on evenly spaced time.
"""
precompute_ridge_slope(args...; kwargs...) = precompute_ridge(args...; kwargs...)[1,:]

"""

    kendalltau(x)

Compute the kendall-τ correlation coefficient of the time series `x`.
"""
kendalltau(x) = corkendall(collect(1.0:length(x)), x)

"""

    spearman(x)

Compute the spearman correlation coefficient of the time series `x`.
"""
spearman(x) = corspearman(collect(1.0:length(x)), x)
