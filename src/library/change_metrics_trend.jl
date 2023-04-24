Base.@kwdef struct ChangeMetricsParams <: Params
    lambda::Real = 0.0
end

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

    RidgeRegression(t, width; lambda = 0.0) → rr

Initialize a ridge regression for a time vector `t`, a sliding-window `width` and
an optional regularizeation term `lambda`. The output `rr` can then be used as a
function to perform the regression and extract the slope!
If `lambda = 0`, linear regression is recovered (default case).
For more information, visit: https://en.wikipedia.org/wiki/Ridge_regression.

# Examples
```jldoctest
julia> t = 0.0:1.0:100;
julia> x = 2.3 .* t .+ 1.2;
julia> rr = RidgeRegression(t);
julia> rr(x)
2.299999999999999
```
"""
struct RidgeRegression{T} <: StructFunction
    equispaced::Bool
    regression_matrix::Matrix{T}
end

function RidgeRegression(t::AbstractVector, cmp::ChangeMetricsParams)
    regression_matrix = precompute_ridge(t, cmp.lambda)
    return RidgeRegression(isequispaced(t), regression_matrix)
end
RidgeRegression(t) = RidgeRegression(isequispaced(t), precompute_ridge(t))

function (rr::RidgeRegression)(x::AbstractVector{T}) where{T<:Real}
    M = rr.regression_matrix
    if !(rr.equispaced)
        error("Time vector is not evenly spaced. So far, the API is only designed for evenly spaced time series!")
        # For future something like: M .= precompute_ridge(window_view(t))
    end
    return view(M, 1, :)' * x     # we are only interested in the slope.
end

function precompute_ridge(t::AbstractVector{T}, lambda::Real) where {T<:Real}
    TT = hcat(t, ones(T, length(t)))'
    return inv(TT * TT' + lambda .* LinearAlgebra.I(2) ) * TT
end