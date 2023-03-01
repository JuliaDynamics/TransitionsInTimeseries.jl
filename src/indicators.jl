"""

    ar1_whitenoise(x)

Estimate the AR1 regression coefficient `θ` of a vector time series `x`.
Computation based on the analytic solution of the least-square parameter estimation:

```math
    \\hat{\\theta} = \\dfrac{\\sum_{i=2}^{n} x_i \\, x_{i-1}}{\\sum_{i=2}^{n} x_{i-1}^2}
```

"""
function ar1_whitenoise(x::AbstractVector{T}) where {T<:Real}
    n = length(x)
    return (view(x, 2:n)' * view(x, 1:n-1)) / (view(x, 1:n-1)' * view(x, 1:n-1))
end

"""

    variance(x)

Compute variance of `x`.
"""
function variance(args...; kwargs...)
    return var(args...; kwargs...)
end

"""

    skw(x)

Compute skewness of `x`.
"""
function skw(args...; kwargs...)
    return skewness(args...; kwargs...)
end

"""

    krt(x)

Compute kurtosis of `x`.
"""
function krt(args...; kwargs...)
    return kurtosis(args...; kwargs...)
end