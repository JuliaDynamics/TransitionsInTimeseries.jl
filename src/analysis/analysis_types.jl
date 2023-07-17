"""
    IndicatorsConfig(t::AbstractVector, indicators::Vector{<:Function}; kwargs...)
    IndicatorsConfig(t::AbstractVector, indicators...; kwargs...)

A configuration for computing indicators from a timeseries with time `t`.
Indicators are standard Julia functions that input an `::AbstractVector` and output
a `::Real`. Any number of indicators can be given, either as a vector or as
variable arguments. Some indicators typically used in
the literature are listed in the documentation section [Indicators](@ref).

this configuration with [`SignificanceConfig`](@ref) this configuration is given to
[`indicators_analysis`](@ref) for the main processing of TransitionsInTimeseries.jl.

The keywords `width, stride` are propagated into [`WindowViewer`](@ref) to create a
sliding window for estimating the indicators.
The time vector corresponding to the indicators timeseries
is obtained from `t` using the keyword `whichtime`. Options include:

- `last`: use the last timepoint of each window
- `midpoint`: use the mid timepoint of each time window
- `first`: use first timepoint of each window

In fact, the indicators time vector is computed simply via
```julia
t_indicator = windowmap(whichtime, t; width, stride)
```
so any other function of the timewindow may be given to extract the time point itself,
such as `mean`.
"""
struct IndicatorsConfig{F}
    t_indicator::AbstractVector
    indicators::Vector{F}
    width::Int
    stride::Int
end

function IndicatorsConfig(
        t::AbstractVector,
        indicators::Vector;
        width = default_window_width(t),
        stride = DEFAULT_WINDOW_STRIDE,
        whichtime =  midpoint,
    )
    indicators = precompute_metrics(indicators, t[1:width])
    t_indicator = windowmap(whichtime, t; width = width, stride = stride)
    return IndicatorsConfig(t_indicator, indicators, width, stride)
end

function IndicatorsConfig(t::AbstractVector, f::Vararg; kwargs...)
    return IndicatorsConfig(t, f_windowmaptime, collect(f), kwargs...)
end

"""
    SignificanceConfig(indconfig, f, change_metrics...; kwargs...)

A configuration for estimating a significant change of indicators based on a sliding
window and and `indconfig::IndicatorsConfig`. The time vector resulting from the
sliding window is obtained under the hood by:
```julia
t_change = windowmap(f, indconfig.t_indicator, kwargs...)
```

Common choices for `f` are:
 - `last`: use `x[k-width:k]` to estimate the indicator change at time step `k`.
 - `midpoint`: use `x[k-width÷2:k+width÷2]` to estimate the indicator change
at time step `k`.
 - `first`: use `x[k:k+width]` to estimate the indicator change at time step `k`.

Along with [`IndicatorsConfig`](@ref), `SignificanceConfig` is given to
[`indicators_analysis`](@ref).

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
- `width::Int`: width given to the [`WindowViewer`](@ref) of the indicator timeseries
  to estimate transient changes of the indicator.
- `stride::Int`: stride given to the [`WindowViewer`](@ref) of the indicator timeseries
  to estimate transient changes of the indicator.
- `tail::Symbol`: kind of tail test to do (one of `:left, :right, :both`).
"""
struct SignificanceConfig{M<:Vector{<:Function}, S<:Surrogate, R<:AbstractRNG}
    t_change::AbstractVector
    change_metrics::M
    n_surrogates::Int
    surrogate_method::S
    rng::R
    width::Int
    stride::Int
    tail::Symbol
end

function SignificanceConfig(
    indconfig::IndicatorsConfig,
    f_windowmaptime::Function,
    change_metrics::Vector;
    n_surrogates::Int = DEFAULT_N_SURROGATES,
    surrogate_method::S = DEFAULT_SURROGATE_METHOD,
    rng::R = Random.default_rng(),
    width::Int = default_window_width(indconfig.t_indicator),
    stride::Int = DEFAULT_WINDOW_STRIDE,
    tail::Symbol = :both,
) where {S<:Surrogate, R<:AbstractRNG}

    change_metrics = precompute_metrics(change_metrics, indconfig.t_indicator[1:width])
    t_change = windowmap(f_windowmaptime, indconfig.t_indicator;
        width = width, stride = stride)

    return SignificanceConfig(t_change, change_metrics, n_surrogates, surrogate_method,
        rng, width, stride, tail)
end

function SignificanceConfig(indconfig::IndicatorsConfig, f_windowmaptime::Function,
    change_metrics::Vararg{Function}; kwargs...)
    return SignificanceConfig(indconfig, f_windowmaptime,
        collect(change_metrics); kwargs...)
end

"""
    IndicatorsResults

A struct containing the output of [`indicators_analysis`](@ref), which is
the main computational part of TransitionsInTimeseries.jl.
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