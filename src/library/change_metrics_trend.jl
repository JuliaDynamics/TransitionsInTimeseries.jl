#####################################################
# Trend metrics
#####################################################
using StatsBase: corkendall, corspearman

"""
    kendalltau(x)

Compute the kendall-τ correlation coefficient of the time series `x`.
"""
kendalltau(x) = corkendall(1:length(x), x)

"""
    spearman(x)

Compute the spearman correlation coefficient of the time series `x`.
"""
spearman(x) = corspearman(1:length(x), x)

"""
    RidgeRegressionSlope(x::AbstractVector)

Returns the slope of the [ridge regression](https://en.wikipedia.org/wiki/Ridge_regression) of `x`.

The computation is based on `lambda_ridge::Real`, a regularizing constant.
"""
Base.@kwdef struct RidgeRegressionSlope <: PrecomputableFunction
    lambda::Real = 0.0
end

struct PrecomputedRidgeRegressionSlope{T} <: Function
    equispaced::Bool
    regression_matrix::Matrix{T}
end

function precompute(rr::RidgeRegressionSlope, t::AbstractVector{T}) where {T<:Real}
    regression_matrix = precompute_ridgematrix(t, rr.lambda)
    return PrecomputedRidgeRegressionSlope(isequispaced(t), regression_matrix)
end

function (rr::PrecomputedRidgeRegressionSlope)(x::AbstractVector{T}) where {T<:Real}
    if !(rr.equispaced)
        error("Time vector is not evenly spaced." * 
            "So far, the API is only designed for evenly spaced time series!")
        # For future something like: M .= precompute_ridgematrix(window_view(t))
    end
    return view(rr.regression_matrix, 1, :)' * x    # only return slope.
    # view(rr.regression_matrix, 2, :)' * x --> would return the bias.
end

function precompute_ridgematrix(t::AbstractVector{T}, lambda::Real) where {T<:Real}
    TT = hcat(t, ones(T, length(t)))'
    return inv(TT * TT' + lambda .* LinearAlgebra.I(2) ) * TT
end