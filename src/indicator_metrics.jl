#####################################################
# Trend metrics
#####################################################
"""

    ridge_regression(x, y)

Perform ridge regression of `y` over `x` with regularization parameter `lambda`.
If `lambda = 0`, linear regression is recovered (default case).
For more information, visit: https://en.wikipedia.org/wiki/Ridge_regression
"""
function ridge_regression(
    x::AbstractVector{T},
    y::AbstractVector{T};
    lambda::T = T(0),
) where {T<:Real}
    X = hcat(x, ones(length(x)))'
    return inv(X * X' + lambda .* I(2)) * X * y
end

#####################################################
# Change-point metrics
#####################################################