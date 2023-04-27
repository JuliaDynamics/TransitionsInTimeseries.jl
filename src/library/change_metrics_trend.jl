#####################################################
# Trend metrics
#####################################################

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
This parameter is passed within a parameter set [`ChangeMetricsParams`](@ref) at
initialization time. It is handled automatically within [`SignificanceConfig`](@ref).
"""
struct RidgeRegressionSlope{T} <: StructFunction
    equispaced::Bool
    regression_matrix::Matrix{T}
end

function RidgeRegressionSlope(t::AbstractVector, cmp::ChangeMetricsParams)
    regression_matrix = precompute_ridge(t, cmp.lambda_ridge)
    return RidgeRegressionSlope(isequispaced(t), regression_matrix)
end
RidgeRegressionSlope(t) = RidgeRegressionSlope(isequispaced(t), precompute_ridge(t))

function (rr::RidgeRegressionSlope)(x::AbstractVector{T}) where{T<:Real}
    M = rr.regression_matrix
    if !(rr.equispaced)
        error("Time vector is not evenly spaced." * 
            "So far, the API is only designed for evenly spaced time series!")
        # For future something like: M .= precompute_ridge(window_view(t))
    end
    return view(M, 1, :)' * x     # we are only interested in the slope.
end

function precompute_ridge(t::AbstractVector{T}, lambda_ridge::Real) where {T<:Real}
    TT = hcat(t, ones(T, length(t)))'
    return inv(TT * TT' + lambda_ridge .* LinearAlgebra.I(2) ) * TT
end

#####################################################
# Change-point metrics
#####################################################