"""
    IndicatorsConfig(indicators...; window_kwargs...)

A configuration for computing indicators from timeseries. Any number of indicators
(standard Julia functions) can be given. Indicators typically used in the literature
are listed in the documentation section [Indicators](@ref).

Keywords are propagated into [`WindowViewer`](@ref)
to create a sliding window for estimating the indicators.

Along with [`SignificanceConfig`](@ref) it is given to [`indicators_analysis`](@ref).
"""
struct IndicatorsConfig{F<:Function, W<:NamedTuple}
    indicators::Vector{<:F}
    window_kwargs::W
end
IndicatorsConfig(f::Vararg{Function}; kwargs...) = IndicatorsConfig(collect(f), kwargs)
IndicatorsConfig(f::Vector{<:Function}; kwargs...) = IndicatorsConfig(f, kwargs)

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
Base.@kwdef struct SignificanceConfig{M<:Vector{<:Function}, S<:Surrogate, R<:AbstractRNG}
    change_metrics::M = [spearman]
    n_surrogates::Int = 10_000
    surrogate_method::S = RandomFourier()
    rng::R = Random.default_rng()
    width::Int = 20
    stride::Int = 1
end

SignificanceConfig(change_metrics::Vararg{Function}; kwargs...) =
SignificanceConfig(collect(change_metrics); kwargs...)

SignificanceConfig(change_metrics::Vector{<:Function}; kwargs...) =
SignificanceConfig(; kwargs..., change_metrics)

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
struct IndicatorsResults{TT, T<:Real, X<:Real, XX<:AbstractVector{X}, F<:Function, Z<:Function}
    t::TT # original time vector; most often it is `Base.OneTo`.
    x::XX

    indicators::Vector{F}
    t_indicator::Vector{T}
    x_indicator::Matrix{X}

    change_metrics::Vector{Z}
    t_change::Vector{T}
    x_change::Matrix{X}
    s_change::Array{X, 3}
end
