# Indicators already in other packages are re-exported
using Statistics: mean, var
using StatsBase: skewness, kurtosis

abstract type Params end
abstract type StructFunction <: Function end
Base.@kwdef struct IndicatorsParams <: Params
    q_lofreq::Real = 0.1
end

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
"""
struct LowfreqPowerSpectrum <: StructFunction
    plan::Any
    i_lofreq::Int
end

function LowfreqPowerSpectrum(x::AbstractVector, iparams::IndicatorsParams)
    p = plan_rfft(x)
    fft_x = p*x
    i_lofreq = Int(round(length(fft_x) * iparams.q_lofreq))
    return LowfreqPowerSpectrum(p, i_lofreq)
end

function (lfps::LowfreqPowerSpectrum)(x::AbstractVector)
    ps = abs.(lfps.plan * x)
    return sum(view(ps, 1:lfps.i_lofreq)) / sum(ps)
end