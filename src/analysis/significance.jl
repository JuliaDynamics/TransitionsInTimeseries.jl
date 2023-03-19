using Statistics: quantile

"""
    SignificanceTest

Abstract type encompassing ways to test for significance in the function
[`significant`](@ref). Using its subtypes, the _real value_
derived from input data
is compared with the _distribution of the surrogate values_.

Subtypes:

- [`Quantile`](@ref)
- [`Sigma`](@ref)
"""
abstract type Significance end

"""
    indicators_significance(res::IndicatorsResults, q::SignificanceTest) → out

Given the output of [`indicators_analysis`](@ref), estimate when the computed
indicators have a statistically significant change, with "significant" defined by
`q`, which is a subtype of [`SignificanceTest`](@ref).

The output is returned as an `n×m` matrix with `n` the length of the time
vector of the change metric timeseries, and `m` the number of change metrics
(which is the same as the number of indicators).
The values of the matrix are Booleans (`true` for significant).

You can use `prod(out; dim = 2)` for logical `&` product across indicators.
"""
function indicators_significance(res::IndicatorsResults, q::Significance)
    out = zeros(Bool, length(res.t_change), size(res.x_change, 2))
    for i in 1:length(res.t_change) # loop over time
        for j in 1:size(res.x_change, 2) # loop over change metrics
            val = res.x_change[i, j] # value of real change metric
            s_vals = view(res.s_change, i, j, :) # values of all surrogates
            out[i, j] = significant(val, s_vals, q)
        end
    end
    return out
end


#####################################################
# Thresholding methods
#####################################################


"""
    significant(v::Real, s_vals::AbstractVector, q::Significance) → true/false

Return `true` if the value `v` is significant when compared to distribution of
surrogate values `s_vals` according to the criterion `q` (see [`SignificanceTest`](@ref)).
If not significant, return `false`.

If given `q::Real`, then [`Quantile`](@ref) is used.
"""
significant(v, s_vals, q::Real) = significant(v, s_vals, Quantile(q))

"""
    Quantile(q = 0.99, dir::Symbol = :updown) <: SignificanceTest

Significance is claimed by comparing the real value `v` with the `q`-th quantile of the
surrogate distribution. Three possibilities are provided depending on `dir`:
- `dir = :updown`: Significance is claimed if `v` is more than the `q`-th or is
  less than the `1-q`-th quantile
- `dir = :up`: if `v` is more than `q`-th
- `dir = :down`: if `v` is less than `1-q`-th
"""
struct Quantile <: Significance
    q::Float64
    dir::Symbol
end
Quantile(q = 0.99, dir = :updown) = Quantile(q, dir)

function significant(val, s_vals, quant::Quantile)
    q = quant.q
    if quant.dir == :updown
        low, high = quantile(s_vals, (1-q, q))
        return val < low || val > high
    elseif quant.dir == :up
        high = quantile(s_vals, q)
        return val > high
    elseif quant.dir == :down
        low = quantile(s_vals, 1-q)
        return val < low
    else
        error()
    end
end




"""
    Sigma(n::Real = 3, dir::Symbol = :updown) <: SignificanceTest

Significance is claimed by comparing the real value `v` with the
mean `μ` and standard deviation `σ` of the surrogate distribution.
Three possibilities are provided depending on `dir`:
- `dir = :updown`: significance is claimed if `v` is more than `μ + n*σ` or less than
  `μ - n*σ`.
- `dir = :up`: if `v` is more than `μ + n*σ`
- `dir = :down`: if `v` is less than `μ - n*σ`
"""
struct Sigma <: Significance
    n::Float64
    dir::Symbol
end
Sigma(n = 3, dir = :updown) = Sigma(n, dir)

using StatsBase: mean_and_std

function significant(val, s_vals, sigma::Sigma)
    n = sigma.n
    μ, σ = mean_and_std(s_vals)
    low, high = μ - n*σ, μ + n*σ
    if sigma.dir == :updown
        return val < low || val > high
    elseif sigma.dir == :up
        return val > high
    elseif sigma.dir == :down
        return val < low
    else
        error()
    end
end

# TODO: everything below is old Jan code
#=
#####################################################
# Thresholding methods
#####################################################

"""

    threshold_indicators(t, significances; n_indicators)

Return the time steps at which `n_indicators` are significant, based on the
`significance` over time `time`.
"""
function threshold_indicators(
    t::Vector{T},
    significances::Array{T, 3};
    n_indicators::Int = size(significances, 3),
) where {T<:AbstractFloat}
    significant_indicators = vec(sum(significances .>= 1, dims = 3))
    return t[significant_indicators .>= n_indicators]
end

#####################################################
# Significance computation
#####################################################

"""

    measure_significance(res, significance_metrics)

Compute some `significance_metrics` for the `IndicatorsResults` output by
[`analyze_indicators`](@ref analyze_indicators).
"""
function measure_significance(
    res::IndicatorsResults{T},
    significance_metrics::Vector{Function},
) where {T<:AbstractFloat}
    significances = similar(res.X_evolution)
    for i in axes(res.X_evolution, 3)
        significances[:, :, i] = measure_significance(
            res.X_evolution[1, :, i],
            res.S_evolution[:, :, i],
            significance_metrics[i],
        )
    end
    return significances
end

function measure_significance(
    res::IndicatorsResults{T},
    significance_metrics::Function,
) where {T<:AbstractFloat}
    sm = repeat(Function[significance_metrics], outer = size(res.X_evolution, 3))
    return measure_significance(res, sm)
end

function measure_significance(
    x_evolution::AbstractVector{T},
    s_evolution::AbstractMatrix{T},
    significance_metric::Function,
) where {T<:AbstractFloat}
    significance = similar(x_evolution)
    for j in eachindex(significance)
        significance[j] = significance_metric(
            x_evolution[j],
            view(s_evolution, :, j),
        )
    end
    return significance
end

#####################################################
# Percentile metrics
#####################################################
"""

    percentile_idx(s; p, symmetric)

If `symmetric == false`, return indices of `0`-th and `p*100`-th percentile of `s`.
If `symmetric == true` and `p=0.95`, return indices of `2.5`-th and `97.5`-th percentile.
"""
function percentile_idx(
    s::AbstractVector{T};
    p::Real = 0.95,
    symmetric::Bool = false,
) where{T<:AbstractFloat}
    n = length(s)
    sorted_idx = sortperm(s)
    i_lo = symmetric ?
        sorted_idx[intround(n * (0.5 - p/2))] :
        sorted_idx[1]
    i_hi = symmetric ?
        sorted_idx[intround(n * (0.5 + p/2))] :
        sorted_idx[intround(n * p)]
    return i_lo, i_hi
end

"""

    normalized_percentile(x, s; p, symmetric)

Compute the normalized percentile of a value `x` within a dataset `s` and w.r.t. a
percentile value `p`.
"""
function normalized_percentile(
    x::T,
    s::AbstractVector{T};
    p::T = T(0.95),
    symmetric::Bool = false,
) where{T<:AbstractFloat}
    i_lo, i_hi = percentile_idx(s, p=p, symmetric=symmetric)
    percentile_center = 0.5 * (s[i_lo] + s[i_hi])
    interpercentile_distance = 0.5 * (s[i_hi] - s[i_lo])
    return (x - percentile_center) / abs(interpercentile_distance)
end

"""

    confidence_interval(x, s; nstd)

Compute wheter an input value `x` is inside (< 1) or outside (> 1) of the confidence
interval of a dataset `s`, defined by `nstd` standard-deviations.
"""
function confidence_interval(
    x::T,
    s::AbstractVector{T};
    nstd::T=2.0,
) where{T<:AbstractFloat}
    return (x - mean(s)) / (nstd * std(s))
end

"""

    which_percentile(x, s)

Compute the percentile represented by the input value `x` w.r.t. the vector `s`.
"""
function which_percentile(
    x::T,
    s::AbstractVector{T},
) where{T<:AbstractFloat}
    return sum(x .> s) / length(s)
end

intround(x::Real) = Int(round(x))
=#