
"""
    LowfreqPowerSpectrum(; q_lofreq = 0.1)

Return a [`PrecomputableFunction`](@ref) containing all the necessary fields to
generate a [`PrecomputedLowfreqPowerSpectrum`](@ref). The latter can be
initialized by [`precompute`](@ref):

```julia
lfps = precompute( LowfreqPowerSpectrum() )
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
    PrecomputedLowfreqPowerSpectrum

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
