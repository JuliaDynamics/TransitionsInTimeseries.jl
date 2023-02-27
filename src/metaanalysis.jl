#=
ns: number of surrogates
nt_indicator: number of data points in the indicator time series
nt_evolution: number of data points in the evolution metric time series
ni: number of indicators to compute
=#
struct IndicatorEvolutionResults{T<:AbstractFloat}
    t_indicator::Vector{T}      # nt_indicator
    X_indicator::Array{T, 3}    # ns=1 x nt_indicator x ni
    t_evolution::Vector{T}      # nt_evolution
    X_evolution::Array{T, 3}    # ns=1 x nt_evolution x ni
    S_evolution::Array{T, 3}    # ns x nt_evolution x ni
end

struct MetaAnalysisParameters{S<:Surrogate, R<:AbstractRNG}
    n_surrogates::Int
    surrogate_method::S
    rng::R
    wv_indicator_width::Int
    wv_indicator_stride::Int
    wv_evolution_width::Int
    wv_evolution_stride::Int
end

"""
    init_metaanalysis_params(kwargs...)

Initialize a `MetaAnalysisParameters` struct with default choices. Custom values can
be set by providing keyword arguments:

- `n_surrogates`: number of surrogates to generate.
- `surrogate_method`: surrogate generation method.
- `rng`: random number generator.
- `wv_indicator_width`: window width for `indicator` computation.
- `wv_indicator_stride`: window stride for `indicator` computation.
- `wv_evolution_width`: window width for `evolution_metric` computation.
- `wv_evolution_stride`: window stride for `evolution_metric` computation.

"""
function init_metaanalysis_params(;
    n_surrogates::Int = 10_000,
    surrogate_method::S = RandomFourier(),
    rng::AbstractRNG = Random.default_rng(),
    wv_indicator_width::Int = 100,
    wv_indicator_stride::Int = 5,
    wv_evolution_width::Int = 50,
    wv_evolution_stride::Int = 5,
) where {S<:Surrogate}
    return MetaAnalysisParameters(
        n_surrogates, surrogate_method, rng,
        wv_indicator_width, wv_indicator_stride,
        wv_evolution_width, wv_evolution_stride,
    )
end
# TODO: init struct based on dimensions of input time-series

"""
    analyze_indicators(t, x, indicators, evolution_metrics, p)

Compute the `indicators` and their `evolution_metrics` for a timeseries `t`, `x` and
its surrogates. If `t` is not provided, it is simply assumed to be `1:length(x)`.
This meta-analysis is performed based on `p::MetaAnalysisParameters`.
"""
function analyze_indicators(
    t::AbstractVector{T},
    x::AbstractVector{T},
    indicators::Vector{Function},
    evolution_metrics::Vector{Function},
    p::MetaAnalysisParameters,
) where {T<:AbstractFloat}

    n_ind = length(indicators)
    indicator_length = get_mapwindowview_length(x,
        p.wv_indicator_width, p.wv_indicator_stride)
    evolution_length = get_mapwindowview_length(1:indicator_length,
        p.wv_evolution_width, p.wv_evolution_stride)

    X_indicator = fill(T(0), 1, indicator_length, n_ind)
    X_evolution = fill(T(0), 1, evolution_length, n_ind)
    t_indicator = fill(T(0), indicator_length)
    t_evolution = fill(T(0), evolution_length)

    for i in 1:n_ind
        if i == 1       # t_indicator and t_evolution only need to be computed once.
            t1, x_indicator, t2, x_evolution = analyze_indicator(
                t, x, indicators[i], evolution_metrics[i], p)
            copy!(t_indicator, t1)
            copy!(t_evolution, t2)
        else
            x_indicator, x_evolution = analyze_indicator(
                x, indicators[i], evolution_metrics[i], p)
        end
        X_indicator[1, :, i] .= x_indicator
        X_evolution[1, :, i] .= x_evolution
    end

    sgen = surrogenerator(x, p.surrogate_method, p.rng)
    S_evolution = fill(T(0.0), p.n_surrogates, evolution_length, n_ind)
    for j in 1:p.n_surrogates
        for i in 1:n_ind
            s = sgen()
            s_evolution = analyze_indicator(
                s, indicators[i], evolution_metrics[i], p)
            S_evolution[j, :, i] .= s_evolution
        end
    end
    return IndicatorEvolutionResults(
        t_indicator, X_indicator, t_evolution, X_evolution, S_evolution,
    )
end

# allow for a single indicator and a single evolution metric.
function analyze_indicators(
    t::AbstractVector{T},
    x::AbstractVector{T},
    indicators::Function,
    evolution_metrics::Function,
    p::MetaAnalysisParameters,
) where {T<:AbstractFloat}
    return analyze_indicators(t, x, Function[indicators], Function[evolution_metrics], p)
end

# allow for multiple indicators and a single evolution metric.
function analyze_indicators(
    t::AbstractVector{T},
    x::AbstractVector{T},
    indicators::Vector{Function},
    evolution_metrics::Function,
    p::MetaAnalysisParameters,
) where {T<:AbstractFloat}
    em = repeat(Function[evolution_metrics], outer = length(indicators))
    return analyze_indicators(t, x, indicators, em, p)
end

# allow the user to not provide any time vector.
function analyze_indicators(
    x::AbstractVector{T},
    indicators::VFi,
    evolution_metrics::VFe,
    p::MetaAnalysisParameters,
) where {
    T<:AbstractFloat,
    F<:Function,
    VFi<:Union{F, Vector{F}},
    VFe<:Union{F, Vector{F}},
}
    t = 1:length(x)
    return analyze_indicators(t, x, indicators, evolution_metrics, p)
end

"""

    analyze_indicator(t, x, indicator, evolution)

Compute a single `indicator` and its `evolution_metric` for a timeseries `t`, `x` and
some surrogates. If `t` is not provided, it is simply assumed to be `1:length(x)`.
This meta-analysis is performed based on `p::MetaAnalysisParameters`.
"""
function analyze_indicator(
    t::AbstractVector{T},
    x::Vector{T},
    indicator::Function,
    evolution_metric::Function,
    p::MetaAnalysisParameters,
) where {T<:AbstractFloat}
    t_indicator, x_indicator = mapwindow(t, x, indicator,
        p.wv_indicator_width, p.wv_indicator_stride)
    t_evolution, x_evolution = mapwindow(t_indicator, x_indicator, evolution_metric,
        p.wv_evolution_width, p.wv_evolution_stride)

    return t_indicator, x_indicator, t_evolution, x_evolution
end

function analyze_indicator(
    x::Vector{T},
    indicator::Function,
    evolution_metric::Function,
    p::MetaAnalysisParameters,
) where {T<:AbstractFloat}
    x_indicator = mapwindow(x, indicator, p.wv_indicator_width, p.wv_indicator_stride)
    x_evolution = mapwindow(x_indicator, evolution_metric,
        p.wv_evolution_width, p.wv_evolution_stride)

    return x_evolution
end

"""

    mapwindow(x, f, wv_width, wv_stride)
    mapwindow(t, x, f, wv_width, wv_stride)

Generate a `WindowViewer` of `x` with `wv_width` and `wv_stride` and map function `f`
over it. If the time vector `t` is provided, additionally return the time vector resulting
from applying the `WindowViewer`.
"""
function mapwindow(
    x::Vector{T},
    f::Function,
    wv_width::Int,
    wv_stride::Int,
) where {T<:AbstractFloat}
    wv = WindowViewer(x, wv_width, wv_stride)
    return map(f, wv)
end

function mapwindow(
    t::AbstractVector{T},
    x::Vector{T},
    f::Function,
    wv_width::Int,
    wv_stride::Int,
) where {T<:AbstractFloat}
    wv = WindowViewer(x, wv_width, wv_stride)
    return t[wv.strided_indices], map(f, wv)
end

"""

    get_mapwindowview_length(x, wv_width, wv_stride)

Compute the length of the `WindowViewer` induced by `x`, `wv_width` and `wv_stride`.
"""
function get_mapwindowview_length(x, wv_width, wv_stride)
    wv = WindowViewer(x, wv_width, wv_stride)
    return length(wv)
end

#####################################################
# Change-point metrics
#####################################################