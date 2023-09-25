using LinearAlgebra: ⋅
using Statistics: var

"""
    ar1_whitenoise(x::AbstractVector)

Return the AR1 regression coefficient of a time series `x` by computing
the analytic solution of the least-square parameter estimation under white-noise
assumption for the data-generating process:

```math
\\theta = \\sum_{i=2}^{n} x_i  x_{i-1} / \\sum_{i=2}^{n} x_{i-1}^2
```
"""
function ar1_whitenoise(x::AbstractVector{<:Real})
    n = length(x)
    return (view(x, 2:n) ⋅ view(x, 1:n-1)) / (view(x, 1:n-1) ⋅ view(x, 1:n-1))
end