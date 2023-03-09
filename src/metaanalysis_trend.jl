#####################################################
# Trend metrics
#####################################################
struct RidgeRegression{T}
    equispaced::Bool
    regression_matrix::Matrix{T}
end

"""

RidgeRegression(t, width; lambda = 0.0) -> rr

Initialize a ridge regression for a time vector `t`, a sliding-window `width` and
an optional regularizeation term `lambda`. The output `rr` can then be used as a
function to perform the regression and extract the slope!

===============
Example
===============
```julia
t = 0.0:1.0:100
m, p = 2.3, 1.2
x = m .* t + p
rr = RidgeRegression(t)
rr(x) ≈ m
```
"""
function RidgeRegression(t::AbstractVector, width::Int; lambda = 0.0)
    equispaced, mean_dt = is_equispaced(t)
    t_regression = range(0.0, step = mean_dt, length = width)
    regression_matrix = precompute_ridge(t_regression, lambda = lambda)
    return RidgeRegression(equispaced, regression_matrix)
end
RidgeRegression(t) = RidgeRegression(is_equispaced(t)[1], precompute_ridge(t))

function (rr::RidgeRegression)(x::AbstractVector{T}) where{T<:Real}
    M = rr.regression_matrix
    if !(rr.equispaced)
        error("Time vector is not evenly spaced. So far, the API is only designed for evenly spaced time series!")
        # For future something like: M .= precompute_ridge(window_view(t))
    end
    return M[1, :]' * x     # we are only interested in the slope.
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

    kendalltau(x)

Compute the kendall-τ correlation coefficient of the time series `x`.
"""
kendalltau(x) = corkendall(collect(1.0:length(x)), x)

"""

    spearman(x)

Compute the spearman correlation coefficient of the time series `x`.
"""
spearman(x) = corspearman(collect(1.0:length(x)), x)
