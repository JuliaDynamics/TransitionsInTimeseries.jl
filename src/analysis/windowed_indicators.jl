"""
    WindowConfig

Supertype used to define the window configuration of [`estimate_indicator_changes`](@ref).
Valid subtypes are:
 - [`SlidingWindowConfig`](@ref).
 - [`SegmentWindowConfig`](@ref).
"""
abstract type WindowConfig end

"""
    SlidingWindowConfig(indicators, change_metrics; kwargs...) → config

A configuration struct for TransitionsInTimeseries.jl that collects the sliding window
parameters, the indicators and the corresponding metrics to use in
[`estimate_indicator_changes`](@ref).

`indicators` is a tuple of indicators (or a single indicator).
`change_metrics` is also a tuple or a single function. If a single function,
the same change metric is used for all provided indicators. This way the analysis
can be efficiently repeated for many indicators and/or change metrics.

Both indicators and change metrics are generic Julia functions that input an
`x::AbstractVector` and output an `s::Real`. Any appropriate function may be given and
see [making custom indicators/change metrics](@ref own_indicator) in the documentation
for more information.

## Keyword arguments

- `width_ind::Int=100, stride_ind::Int=1`: width and stride given to [`WindowViewer`](@ref)
  to compute the indicator from the input timeseries.
- `width_cha::Int=50, stride_cha::Int=1`: width and stride given to [`WindowViewer`](@ref)
  to compute the change metric timeseries from the indicator timeseries.
- `whichtime = midpoint`: The time vector corresponding to the indicators / change metric
  timeseries is obtained from `t` in [`estimate_indicator_changes`](@ref) using the keyword
  `whichtime`. Options include:
    - `last`: use the last timepoint of each window
    - `midpoint`: use the mid timepoint of each time window
    - `first`: use first timepoint of each window
  In fact, the indicators time vector is computed simply via
  ```julia
  t_indicator = windowmap(whichtime, t; width_ind, stride_ind)
  t_change = windowmap(whichtime, t_indicator; width_cha, stride_cha)
  ```
  so any other function of the time window may be given to extract the time point itself,
  such as `mean` or `median`.

- `T = Float64`: Element type of input timeseries to initialize some computations.

"""
struct SlidingWindowConfig{F, G, W<:Function} <: WindowConfig
    indicators::F
    change_metrics::G
    width_ind::Int
    stride_ind::Int
    width_cha::Int
    stride_cha::Int
    whichtime::W
end

function SlidingWindowConfig(
        indicators, change_metrics;
        width_ind = 100,
        stride_ind = 1,
        width_cha = 50,
        stride_cha = 1,
        whichtime =  midpoint,
        T = Float64,
    )
    indicators, change_metrics = sanitycheck_metrics(indicators, change_metrics)
    # Last step: precomputable functions, if any
    indicators = map(f -> precompute(f, 1:T(width_ind)), indicators)
    change_metrics = map(f -> precompute(f, 1:T(width_cha)), change_metrics)

    return SlidingWindowConfig(
        indicators, change_metrics,
        width_ind, stride_ind, width_cha, stride_cha, whichtime,
    )
end

"""
    SegmentWindowConfig(indicators, change_metrics, tseg_start, tseg_end;
        kwargs...) → config

A configuration struct for TransitionsInTimeseries.jl that collects the start and end time
of the segments over which the analysis [`estimate_indicator_changes`](@ref) is to be run,
as well as its windowing parameters.

`indicators` is a tuple of indicators (or a single indicator).
`change_metrics` is also a tuple or a single function. If a single function,
the same change metric is used for all provided indicators. This way the analysis
can be efficiently repeated for many indicators and/or change metrics.
`tseg_start` and `tseg_end` are `::Vector` of length `n`, with `n` the number of segments.

Both indicators and change metrics are generic Julia functions that input an
`x::AbstractVector` and output an `s::Real`. Any appropriate function may be given and
see [making custom indicators/change metrics](@ref own_indicator) in the documentation
for more information.

## Keyword arguments

- `width_ind::Int=100, stride_ind::Int=1`: width and stride given to [`WindowViewer`](@ref)
  to compute the indicator from the input timeseries.
- `min_width_cha::Int=50`: minimal width required to perform the change metric estimation.
  If segment not sufficiently long, return `Inf`.
- `whichtime = midpoint`: The time vector corresponding to the indicators / change metric
  timeseries is obtained from `t` in [`estimate_indicator_changes`](@ref) using the keyword
  `whichtime`. Options include:
    - `last`: use the last timepoint of each window
    - `midpoint`: use the mid timepoint of each time window
    - `first`: use first timepoint of each window
  In fact, the indicators time vector for the k-th segment is computed simply via:
  ```julia
  t_indicator[k] = windowmap(whichtime, t; width_ind, stride_ind)
  ```
  so any other function of the time window may be given to extract the time point itself,
  such as `mean` or `median`.

- `T = Float64`: Element type of input timeseries to initialize some computations.

"""
struct SegmentWindowConfig{F, G, W<:Function} <: WindowConfig
    indicators::F
    change_metrics::G
    tseg_start::Vector
    tseg_end::Vector
    width_ind::Int
    stride_ind::Int
    min_width_cha::Int
    whichtime::W
end

function SegmentWindowConfig(
    indicators, change_metrics, tseg_start, tseg_end;
    width_ind = 100,
    stride_ind = 1,
    min_width_cha = 50,
    whichtime =  midpoint,
    T = Float64,
)
    if length(tseg_start) ≠ length(tseg_end)
        throw(ArgumentError("The vectors containing the start and end time of the"*
            " segments must be of equal length."))
    end
    indicators, change_metrics = sanitycheck_metrics(indicators, change_metrics)
    # Last step: precomputable functions, if any
    indicators = map(f -> precompute(f, 1:T(width_ind)), indicators)

    return SegmentWindowConfig(
        indicators, change_metrics, tseg_start, tseg_end,
        width_ind, stride_ind, min_width_cha, whichtime,
    )
end


function sanitycheck_metrics(indicators, change_metrics)
    if !(indicators isa Tuple)
        indicators = (indicators,)
    end
    if !(change_metrics isa Tuple)
        change_metrics = (change_metrics,)
    end
    L = length(indicators)
    if length(change_metrics) ∉ (1, L)
        throw(ArgumentError("The amount of change metrics must be as many as the"*
            "indicators, or only 1."))
    end
    return indicators, change_metrics
end

"""
    estimate_indicator_changes(config::SlidingWindowConfig, x [,t]) → output

Estimate possible transitions for input timeseries `x` using a sliding window approach
as described by `config`:

1. Estimate the timeseries of an indicator by sliding a window over the input timeseries.
2. Estimate changes of an indicator by sliding a window of the change metric over
   the indicator timeseries.

If `t` (the time vector of `x`), is not provided, it is assumed `t = eachindex(x)`.

Return the output as [`WindowResults`](@ref) which can be given to
[`significant_transitions`](@ref) to deduce which possible transitions are statistically
significant using a variety of significance tests.
"""
function estimate_indicator_changes(config::SlidingWindowConfig, x, t = eachindex(x))
    # initialize time vectors
    t_indicator = windowmap(config.whichtime, t; width = config.width_ind,
        stride = config.stride_ind)
    t_change = windowmap(config.whichtime, t_indicator; width = config.width_cha,
        stride = config.stride_cha)
    len_ind = length(t_indicator)
    len_change = length(t_change)

    # initialize array containers
    (; indicators, change_metrics) = config
    X = eltype(x)
    n_ind = length(indicators)
    x_indicator = zeros(X, len_ind, n_ind)
    x_change = zeros(X, len_change, n_ind)
    one2one = length(change_metrics) == length(indicators)

    # Loop over indicators
    for i in 1:n_ind
        indicator = config.indicators[i]
        chametric = one2one ? change_metrics[i] : change_metrics[1]
        z = view(x_indicator, :, i)
        windowmap!(indicator, z, x;
            width = config.width_ind, stride = config.stride_ind
        )
        windowmap!(chametric, view(x_change, :, i), z;
            width = config.width_cha, stride = config.stride_cha
        )
    end

    # put everything together in the output type
    return SlidingWindowResults(
        t, x, t_indicator, x_indicator, t_change, x_change, config
    )
end

"""
    estimate_indicator_changes(config::SegmentWindowConfig, x [,t]) → output

Estimate possible transitions for input timeseries `x` using a segmented window approach
as described by `config`:

1. For each segment specified in `config`, estimate the corresponding indicator timeseries
   by sliding a window over the input timeseries.
2. For each segment of the indicator timeseries, estimate a scalar change metric by applying
   a window of the size of the indicator segment.

If `t` (the time vector of `x`), is not provided, it is assumed `t = eachindex(x)`.

Return the output as [`WindowResults`](@ref) which can be given to
[`significant_transitions`](@ref) to deduce which possible transitions are statistically
significant using a variety of significance tests.
"""
function estimate_indicator_changes(config::SegmentWindowConfig, x, t)
    X, T = eltype(x), eltype(t)
    (; indicators, change_metrics, tseg_start, tseg_end) = config
    n_ind = length(indicators)

    t_indicator = [T[] for _ in eachindex(tseg_start)]
    x_indicator = [X[;;] for _ in eachindex(tseg_start)]
    t_change = T.(config.tseg_end)
    x_change = fill(Inf, length(tseg_start), n_ind)
    one2one = length(change_metrics) == length(indicators)

    for k in eachindex(tseg_start)
        tseg, xseg = segment(t, x, tseg_start[k], tseg_end[k])
        t_indicator[k] = windowmap(config.whichtime, tseg; width = config.width_ind,
            stride = config.stride_ind)
        len_ind = length(t_indicator[k])

        # Init with Inf instead of 0 to easily recognise when the segment was too short
        # for the computation to be performed.
        x_indicator[k] = fill(Inf, len_ind, n_ind)  

        # only analyze if segment long enough to compute metrics
        if len_ind > config.min_width_cha
            # Loop over indicators
            for i in 1:n_ind
                indicator = indicators[i]
                chametric = one2one ? change_metrics[i] : change_metrics[1]
                z = view(x_indicator[k], :, i)
                windowmap!(indicator, z, xseg;
                    width = config.width_ind, stride = config.stride_ind
                )
                if chametric isa PrecomputableFunction
                    chametric = precompute(chametric, t_indicator[k])
                end
                x_change[k, i] = chametric(z)
            end
        end
    end
    # put everything together in the output type
    return SegmentWindowResults(t, x, t_indicator, x_indicator, t_change, x_change, config)
end

function segment(t, x, t1, t2)
    i1, i2 = argmin(abs.(t .- t1)), argmin(abs.(t .- t2))
    return t[i1:i2], x[i1:i2]
end

"""
    WindowResults

Supertype used to gather results of [`estimate_indicator_changes`](@ref).
Valid subtypes are:
 - [`SlidingWindowResults`](@ref).
 - [`SegmentWindowResults`](@ref).
"""
abstract type WindowResults end

"""
    SlidingWindowResults

A struct containing the output of [`estimate_indicator_changes`](@ref) used with
[`SlidingWindowConfig`](@ref). It can be used for further analysis, visualization,
or given to [`significant_transitions`](@ref).

It has the following fields that the user may access

- `x`: the input timeseries.
- `t`: the time vector of the input timeseries.

- `x_indicator`, the indicator timeseries (matrix with each column one indicator).
- `t_indicator`, the time vector of the indicator timeseries.

- `x_change`, the change metric timeseries (matrix with each column one change metric).
- `t_change`, the time vector of the change metric timeseries.

- [`config::SlidingWindowConfig`](@ref), used for the analysis.
"""
struct SlidingWindowResults{TT, T<:Real, X<:Real, XX<:AbstractVector{X},
    W} <: WindowResults
    t::TT # original time vector; most often it is `Base.OneTo`.
    x::XX
    t_indicator::Vector{T}
    x_indicator::Matrix{X}
    t_change::Vector{T}
    x_change::Matrix{X}
    config::W
end

"""
    SegmentWindowResults

A struct containing the output of [`estimate_indicator_changes`](@ref) used with
[`SegmentWindowConfig`](@ref). It can be used for further analysis, visualization,
or given to [`significant_transitions`](@ref).

It has the following fields that the user may access

- `x`: the input timeseries.
- `t`: the time vector of the input timeseries.

- `x_indicator::Vector{Matrix}`, with `x_indicator[k]` the indicator timeseries (matrix
   with each column one indicator) of the `k`-th segment.
- `t_indicator::Vector{Vector}`, with `t_indicator[k]` the time vector of the indicator
  timeseries for the `k`-th segment.

- `x_change::Matrix`, the change metric values with `x[k, i]` the change metric of the
  `i`-th indicator for the `k`-th segment.
- `t_change`, the time vector of the change metric.

- [`config::SegmentWindowConfig`](@ref), used for the analysis.
"""
struct SegmentWindowResults{TT, T<:Real, X<:Real, XX<:AbstractVector{X},
    W} <: WindowResults
    t::TT # original time vector; most often it is `Base.OneTo`.
    x::XX
    t_indicator::Vector{Vector{T}}
    x_indicator::Vector{Matrix{X}}
    t_change::Vector{T}
    x_change::Matrix{X}
    config::W
end

function Base.show(io::IO, ::MIME"text/plain", res::WindowResults)
    println(io, "WindowResults")
    descriptors = [
        "input timeseries" => summary(res.x),
        "indicators" => [nameof(i) for i in res.config.indicators],
        "indicator (window, stride)" => (res.config.width_ind, res.config.stride_ind),
        "change metrics" => [nameof(c) for c in res.config.change_metrics],
        show_changemetric(res),
    ]
    padlen = maximum(length(d[1]) for d in descriptors) + 2
    for (desc, val) in descriptors
        println(io, rpad(" $(desc): ", padlen), val)
    end
end

function show_changemetric(res::SlidingWindowResults)
    return "change metric (window, stride)" => (res.config.width_cha, res.config.stride_cha)
end

function show_changemetric(res::SegmentWindowResults)
    return "change metric (window)" => ([length(t) for t in res.t_indicator])
end