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
    LowfreqPowerSpectrum(x::AbstractVector)

Returns a number between `0` and `1` that characterizes the amount of power contained
in the low frequencies of the power density spectrum of `x`.
A result of `0.1` means that 10% of the integrated power is contained in low frequencies.

The computation is based on `q_lofreq::Real`, a number between `0` and `1` that
determines which frequencies are considered to be "low". This parameter is passed
within a parameter set [`IndicatorsParams`](@ref) at initialization time.
It is handled automatically within [`IndicatorsConfig`](@ref).
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

"""
    PermutationEntropy(x::AbstractVector)

Returns the permutation entropy of `x`. This computation breaks down in computing
the [`entropy_normalized`](https://juliadynamics.github.io/ComplexityMeasures.jl/stable/entropies/#ComplexityMeasures.entropy_normalized) of a [`SymbolicPermutation`](https://juliadynamics.github.io/ComplexityMeasures.jl/stable/probabilities/#ComplexityMeasures.SymbolicPermutation).

The latter is based on `m_perm::Int`, the order of the permutation.
This parameter is passed within a parameter set [`IndicatorsParams`](@ref) at
initialization time. It is handled automatically within [`IndicatorsConfig`](@ref).
"""
struct PermutationEntropy <: StructFunction
    probest::SymbolicPermutation
end

function PermutationEntropy(x::AbstractVector, iparams::IndicatorsParams)
    probest = SymbolicPermutation(m = iparams.m_perm)
    return PermutationEntropy(probest)
end

function (perment::PermutationEntropy)(x::AbstractVector)
    return entropy_normalized(perment.probest, x)
end