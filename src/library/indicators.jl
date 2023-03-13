# Indicators already in other packages are re-exported
using Statistics: mean, var
using StatsBase: skewness, kurtosis

"""
    ar1_whitenoise(x) → ̂θ

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
    midpoint(x)

Return `x[midindex]` with `midindex = round(Int, 0.5(firstindex(x) + lastindex(x)))`.

Typically useful in [`windowmap`](@ref) with a time vector.
"""
midpoint(x) = x[round(Int, 0.5(firstindex(x) + lastindex(x)))]
midpoint(x::Vector) = x[length(x)÷2] # faster