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
    LowfreqPowerSpectrum()

Return a [`PrecomputableFunction`](@ref) containing all the necessary fields to
generate a [`PrecomputedLowfreqPowerSpectrum`](@ref). The latter can be
initialized by [`precompute`](@ref):

```julia
lfps = precompute( LowfreqPowerSpectrum(q_lofreq = 0.1) )
```

Keyword arguments:
 - `q_lofreq`: a number between `0` and `1` that characterises which portion of the
 frequency spectrum is considered to be low. For instance, `q_lofreq = 0.1` implies
 that the lowest 10% of frequencies are considered to be the low ones.
"""
Base.@kwdef struct LowfreqPowerSpectrum <: PrecomputableFunction
    q_lofreq::Real = 0.1
end

"""
    PrecomputedLowfreqPowerSpectrum(x::AbstractVector)

A struct containing all the precomputed fields to efficiently perform repetitive
computation of the low-frequency power spectrum (LFPS), a number between `0` and `1`
that characterizes the amount of power contained in the low frequencies of the power
density spectrum of `x`. Once `lfps::PrecomputedLowfreqPowerSpectrum` is initialized,
it can be used as a function to obtain the LFPS of `x::AbstractVector` by:

```julia
lfps(x)
```
"""
struct PrecomputedLowfreqPowerSpectrum <: Function
    plan::FFTW.cFFTWPlan{ComplexF64, -1, false, 1, UnitRange{Int64}}
    i_pos::Int              # fftfreq(N)[1:i_pos] only includes positive freqs.
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
Base.@kwdef struct PermutationEntropy{S<:SymbolicPermutation} <: Function
    probest::S = SymbolicPermutation(m = 3)
end
PermutationEntropy(m::Int) = PermutationEntropy( SymbolicPermutation(m = m) )

function (perment::PermutationEntropy)(x::AbstractVector)
    return entropy_normalized(perment.probest, x)
end