"""
    IndicatorsConfig(t, indicators...; kwargs...)

A configuration for computing indicators from timeseries with time vector `t`.
Indicators are standard Julia functions that input an `AbstractVector` and output
a real number. Any number of indicators can be given. Indicators typically used in
the literature are listed in the documentation section [Indicators](@ref).

Keywords are propagated into [`WindowViewer`](@ref)
to create a sliding window for estimating the indicators.

Along with [`SignificanceConfig`](@ref) it is given to [`indicators_analysis`](@ref).

## Keyword arguments
- `indicators_params::IndicatorsParams`: contains the parameters related to the
  indicator functions and needed for initialization in [`init_metrics`](@ref)
- `width::Int`: width given to the [`WindowViewer`](@ref) of the input data
  to estimate indicators.
- `stride::Int`: stride given to the [`WindowViewer`](@ref) of the input data
  to estimate indicators.
"""
struct IndicatorsConfig{F<:Function}
    t_indicator::AbstractVector
    indicators::Vector{F}
    indicators_params::IndicatorsParams
    width::Int
    stride::Int
end

function IndicatorsConfig(
    t::AbstractVector,
    indicators::Vector;
    iparams = IndicatorsParams(),
    width = default_window_width(t),
    stride = default_window_stride(t),
    bracketing::Symbol = :left,
)
    indicators = init_metrics(indicators, t[1:width], iparams)
    t_indicator = slidebracket(t, bracketing, width = width, stride = stride)
    return IndicatorsConfig(t_indicator, indicators, iparams, width, stride)
end

function IndicatorsConfig(
    t::AbstractVector,
    f::Vararg{Function};
    kwargs...)
    return IndicatorsConfig(t, collect(f), kwargs...)
end

"""
    SignificanceConfig(indconfig, change_metrics...; kwargs...)

A configuration for estimating a significant change of indicators computed from timeseries
based on `indconfig::IndicatorsConfig`.
Along with [`IndicatorsConfig`](@ref) it is given to [`indicators_analysis`](@ref).

`change_metrics` can be a single Julia function, in which case
the same metric is applied over all indicators in [`IndicatorsConfig`](@ref).
`change_metrics` can also be many functions, each one for each indicator.

## Keyword arguments
- `change_metrics_params::ChangeMetricsParams`
- `n_surrogates::Int = 10_000`: how many surrogates to create.
- `surrogate_method::S = RandomFourier()`: what method to use to create the surrogates.
  Any `Surrogate` subtype from [TimeseriesSurrogates.jl](
    https://juliadynamics.github.io/TimeseriesSurrogates.jl/stable/#Surrogate-methods
  ) is valid.
- `rng::AbstractRNG = Random`.default_rng()`: a random number generator for the surrogates.
  See the [Julia manual](https://docs.julialang.org/en/v1/stdlib/Random/#Random-Numbers) for more.
- `width::Int`: width given to the [`WindowViewer`](@ref) of the indicator timeseries
  to estimate transient changes of the indicator.
- `stride::Int`: stride given to the [`WindowViewer`](@ref) of the indicator timeseries
  to estimate transient changes of the indicator.
- `tail::Symbol`: kind of tail test to do (one of `:left, :right, :both`).
"""
struct SignificanceConfig{M<:Vector{<:Function}, S<:Surrogate, R<:AbstractRNG}
    t_change::AbstractVector
    change_metrics::M
    change_metrics_params::ChangeMetricsParams
    n_surrogates::Int
    surrogate_method::S
    rng::R
    width::Int
    stride::Int
    bracketing::Symbol
    tail::Symbol
end

function SignificanceConfig(
    indconfig::IndicatorsConfig,
    change_metrics::Vector;
    change_metrics_params::ChangeMetricsParams = ChangeMetricsParams(),
    n_surrogates::Int = default_n_surrogates(),
    surrogate_method::S = default_surrogate_method(),
    rng::R = Random.default_rng(),
    width::Int = default_window_width(indconfig.t_indicator),
    stride::Int = default_window_stride(indconfig.t_indicator),
    bracketing::Symbol = :left,
    tail::Symbol = :right,
) where {S<:Surrogate, R<:AbstractRNG}

    change_metrics = init_metrics(change_metrics, indconfig.t_indicator[1:width],
        change_metrics_params)
    t_change = slidebracket(indconfig.t_indicator, bracketing, width = width, stride = stride)

    return SignificanceConfig(
        t_change, change_metrics, change_metrics_params,
        n_surrogates, surrogate_method, rng,
        width, stride, bracketing, tail)
end

function SignificanceConfig(
    indconfig::IndicatorsConfig,
    change_metrics::Vararg{Function};
    kwargs...)
    return SignificanceConfig(indconfig, collect(change_metrics); kwargs...)
end

"""
    IndicatorsResults

A struct containing the output of [`indicators_analysis`](@ref), which is
the main computational part of TransitionIndicators.jl.
It can be used for analysis and visualization purposes.

It has the following fields:

- `x`: the input timeseries.
- `t`: the time vector of the input timeseries.

- `indicators::Vector{Function}`: indicators used in the processing.
- `x_indicator`, the indicator timeseries (matrix with each column one indicator).
- `t_indicator`, the time vector of the indicator timeseries.

- `change_metrics::Vector{Function}`: change metrics used in the processing.
- `x_change`, the change metric timeseries (matrix with each column one change metric).
- `t_change`, the time vector of the change metric timeseries.
- `pval`, the p-value of the change metrics w.r.t. the surrogates.
  It is a 2-dimensional array, where first dimension = time, second dimension = change
  metric. I.e. `pval[:, k]` will give the time series of p-value of the `k`-th change metric.
"""
struct IndicatorsResults{TT, T<:Real, X<:Real, XX<:AbstractVector{X},
        F<:Function, Z<:Function, S<:Surrogate}
    t::TT # original time vector; most often it is `Base.OneTo`.
    x::XX

    indicators::Vector{F}
    t_indicator::Vector{T}
    x_indicator::Matrix{X}

    change_metrics::Vector{Z}
    t_change::Vector{T}
    x_change::Matrix{X}
    pval::Matrix{X}
    surrogate_method::S
end

function Base.show(io::IO, ::MIME"text/plain", res::IndicatorsResults)
    println(io, "IndicatorsResults")
    descriptors = [
        "input timeseries" => summary(res.x),
        "indicators" => res.indicators,
        "indicator length" => length(res.t_indicator),
        "change metrics" => res.change_metrics,
        "surrogate" => res.surrogate_method,
        "p-value" => summary(res.pval),
    ]
    padlen = maximum(length(d[1]) for d in descriptors) + 3
    for (desc, val) in descriptors
        println(io, rpad(" $(desc): ", padlen), val)
    end
end
