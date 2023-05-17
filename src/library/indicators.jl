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
determines which frequencies are considered to be "low".
"""
Base.@kwdef struct LowfreqPowerSpectrum <: PrecomputableFunction
    q_lofreq::Real = 0.1
end
struct PrecomputedLowfreqPowerSpectrum <: Function
    plan::FFTW.FFTWPlan     # plan for FFT
    i_pos                   # fftfreq(N)[1:i_pos] only includes positive freqs.
    i_lofreq::Int           # fftfreq(N)[1:i_lofreq] only includes low freqs (>0).
end

function precompute(lfps::LowfreqPowerSpectrum, t::AbstractVector)
    p = plan_fft(t)
    i_pos = Int(ceil(length(t) / 2))                # c.f. struct LowfreqPowerSpectrum
    i_lofreq = Int(round(i_pos * lfps.q_lofreq))    # c.f. struct LowfreqPowerSpectrum
    return PrecomputedLowfreqPowerSpectrum(p, i_pos, i_lofreq)
end

function (lfps::PrecomputedLowfreqPowerSpectrum)(x::AbstractVector)
    ps = abs.(lfps.plan * x)
    return sum(view(ps, 1:lfps.i_lofreq)) / sum(view(ps, 1:lfps.i_pos))
end

"""
    PermutationEntropy(x::AbstractVector)

Returns the permutation entropy of `x`. This computation breaks down in computing
the [`entropy_normalized`](https://juliadynamics.github.io/ComplexityMeasures.jl/stable/entropies/#ComplexityMeasures.entropy_normalized) of a [`SymbolicPermutation`](https://juliadynamics.github.io/ComplexityMeasures.jl/stable/probabilities/#ComplexityMeasures.SymbolicPermutation).
"""
Base.@kwdef struct PermutationEntropy <: Function
    probest::SymbolicPermutation = SymbolicPermutation(m = 3)
end
PermutationEntropy(m::Int) = PermutationEntropy( SymbolicPermutation(m = m) )

function (perment::PermutationEntropy)(x::AbstractVector)
    return entropy_normalized(perment.probest, x)
end