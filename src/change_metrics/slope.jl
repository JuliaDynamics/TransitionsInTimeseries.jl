#####################################################
# Trend metrics
#####################################################
using StatsBase: corkendall, corspearman

"""
    kendalltau(x::AbstractVector)

Compute the kendall-τ correlation coefficient of the time series `x`.
`kendalltau` can be used as a change metric focused on trend.
"""
kendalltau(x) = corkendall(1:length(x), x)

"""
    spearman(x::AbstractVector)

Compute the spearman correlation coefficient of the time series `x`.
`spearman` can be used as a change metric focused on trend.
"""
spearman(x) = corspearman(1:length(x), x)

"""
    RidgeRegressionSlope(; lambda = 0.0) → rr

Return a [`PrecomputableFunction`](@ref) containing all the necessary fields to
generate a [`PrecomputedRidgeRegressionSlope`](@ref).
`rr` can be used as a change metric focused on trend.

`lambda` is a regularization constant, usually between `0` and `1`.
"""
Base.@kwdef struct RidgeRegressionSlope <: PrecomputableFunction
    lambda::Real = 0.0
end

struct PrecomputedRidgeRegressionSlope{T} <: Function
    equispaced::Bool
    regression_matrix::Matrix{T}
end

function precompute(rr::RidgeRegressionSlope, t::AbstractVector{T}) where {T<:Real}
    regression_matrix = ridgematrix(t, rr.lambda)
    return PrecomputedRidgeRegressionSlope(isequispaced(t), regression_matrix)
end

"""
    PrecomputedRidgeRegressionSlope(x::AbstractVector)

Return the slope of the [ridge regression](https://en.wikipedia.org/wiki/Ridge_regression) of `x`.
"""
function (rr::PrecomputedRidgeRegressionSlope)(x::AbstractVector{<:Real})
    if !(rr.equispaced)
        error("Time vector is not evenly spaced." *
            "So far, the API is only designed for evenly spaced time series!")
        # For future something like: M .= ridgematrix(window_view(t))
    end
    return view(rr.regression_matrix, 1, :)' * x    # only return slope.
    # view(rr.regression_matrix, 2, :)' * x --> would return the bias.
end

function ridgematrix(t::AbstractVector{T}, lambda::Real) where {T<:Real}
    TT = hcat(t, ones(T, length(t)))'
    return inv(TT * TT' + lambda .* LinearAlgebra.I(2) ) * TT
end