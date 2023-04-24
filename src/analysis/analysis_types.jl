function init_structfuntions(
    change_metrics::Vector,
    t::AbstractVector,
    p::Params,
)
    f = Function[]
    for change_metric in change_metrics
        if change_metric in subtypes(StructFunction)
            push!(f, change_metric(t, p))
        elseif isa(change_metric, Function)
            push!(f, change_metric)
        else
            error("The function $change_metric has an erroneous type.")
        end
    end
    return f
end

"""
    IndicatorsConfig(indicators...; window_kwargs...)

A configuration for computing indicators from timeseries.
Indicators are standard Julia functions that input an `AbstractVector` and output
a real number. Any number of indicators can be given.
Indicators typically used in the literature
are listed in the documentation section [Indicators](@ref).

Keywords are propagated into [`WindowViewer`](@ref)
to create a sliding window for estimating the indicators.

Along with [`SignificanceConfig`](@ref) it is given to [`indicators_analysis`](@ref).
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
)
    indicators = init_structfuntions(indicators, t[1:width], iparams)
    return IndicatorsConfig(t[width:stride:end], indicators, iparams, width, stride)
end

function IndicatorsConfig(
    t::AbstractVector,
    f::Vararg{Function};
    kwargs...)
    return IndicatorsConfig(t, collect(f), kwargs...)
end

"""
    SignificanceConfig(change_metrics...; kwargs...)

A configuration for estimating a significant change of indicators computed from timeseries.
Along with [`IndicatorsConfig`](@ref) it is given to [`indicators_analysis`](@ref).

`change_metrics` can be a single Julia function, in which case
the same metric is applied over all indicators in [`IndicatorsConfig`](@ref).
`change_metrics` can also be many functions, each one for each indicator.

## Keyword arguments
- `n_surrogates::Int = 10_000`: how many surrogates to create.
- `surrogate_method::S = RandomFourier()`: what method to use to create the surrogates.
  Any `Surrogate` subtype from [TimeseriesSurrogates.jl](
    https://juliadynamics.github.io/TimeseriesSurrogates.jl/stable/#Surrogate-methods
  ) is valid.
- `rng::AbstractRNG = Random`.default_rng()`: a random number generator for the surrogates.
  See the [Julia manual](https://docs.julialang.org/en/v1/stdlib/Random/#Random-Numbers) for more.
- `wv_indicator_width::Int = 100`,
- `wv_indicator_stride::Int = 5`,
- `width, stride = 20, 1`: width and stride given to the [`WindowViewer`](@ref) of the
  indicator timeseries. Notice that here the default values are different.
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
) where {S<:Surrogate, R<:AbstractRNG}

    change_metrics = init_structfuntions(
        change_metrics,
        indconfig.t_indicator[1:width],
        change_metrics_params)
    return SignificanceConfig(
        indconfig.t_indicator[width:stride:end],
        change_metrics,
        change_metrics_params,
        n_surrogates,
        surrogate_method,
        rng,
        width,
        stride,
    )
end

function SignificanceConfig(
    indconfig::IndicatorsConfig,
    change_metrics::Vararg{Function};
    kwargs...)
    return SignificanceConfig(indconfig, collect(change_metrics); kwargs...)
end

default_n_surrogates() = 10_000
default_surrogate_method() = RandomFourier()

"""
    IndicatorsResults

A struct containing the output of [`indicators_analysis`](@ref), which is
the main computational part of TransitionIndicators.jl.
It can be given to the [`indicators_significance`](@ref) function
or used for visualizations.

It has the following fields:

- `x`: the input timeseries
- `t`: the time vector of the input timeseries

- `indicators::Vector{Function}`: indicators used in the processing
- `x_indicator`, the indicator timeseries (matrix with each column one indicator)
- `t_indicator`, the indicator timeseries time vector

- `change_metrics::Vector{Function}`: change metrics used in the processing
- `x_change`, the change metric timeseries (matrix with each column one change metric)
- `t_change`, the change metric timeseries time vector
- `s_change`, the result of computing the change metrics for the surrogates.
  It is a 3-dimensional array, where first dimension = time, second dimension = change
  metric, and third dimension = surrogate number. I.e.,
  `S_change[:, j, k]` will give the `k`-th surrogate timeseries of the `j`-th change metric.
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
    s_change::Array{X, 3}
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
        "surrogate #" => size(res.s_change, 3),
    ]
    padlen = maximum(length(d[1]) for d in descriptors) + 3
    for (desc, val) in descriptors
        println(io, rpad(" $(desc): ", padlen), val)
    end
end
