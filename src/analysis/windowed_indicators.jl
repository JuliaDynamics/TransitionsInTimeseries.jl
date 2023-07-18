"""
    WindowedIndicatorMetrics(indicators, change_metrics; kwargs...) → config

A configuration struct for TransitionsInTimeseries.jl that collects
what indicators and corresponding metrics to use in the [`transitions_analysis`](@ref).

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
  timeseries is obtained from `t` in [`transitions_analysis`](@ref) using the keyword
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
struct WindowedIndicatorMetrics{F, G, W<:Function}
    indicators::F
    change_metrics::G
    width_ind::Int
    stride_ind::Int
    width_cha::Int
    stride_cha::Int
    whichtime::W
end

function WindowedIndicatorMetrics(
        indicators, change_metrics;
        width_ind = 100,
        stride_ind = 1,
        width_cha = 50,
        stride_cha = 1,
        whichtime =  midpoint,
        T = Float64,
    )
    # Sanity checks
    if !(indicators isa Tuple)
        indicators = (indicators,)
    end
    if !(change_metrics isa Tuple)
        change_metrics = (change_metrics,)
    end
    L = length(indicators)
    if length(change_metrics) ∉ ≠ (1, L)
        throw(ArgumentError("The amount of change metrics must be as many as the indicators, or only 1."))
    end
    # Last step: precomputable functions, if any
    indicators = map(f -> precompute(f, 1:T(width_ind)), indicators)
    change_metrics = map(f -> precompute(f, 1:T(width_cha)), change_metrics)

    return WindowedIndicatorMetrics(
        indicators, change_metrics,
        width_ind, stride_ind, width_cha, stride_cha, whichtime,
    )
end

"""
    estimate_transitions(config::WindowedIndicatorMetrics, x [,t]) → output

Estimate possible transitions for input timeseries `x` using a sliding window approach
as described by `config`:

1. Estimate the timeseries of an indicator by sliding a window over the input timeseries.
2. Estimate changes of an indicator by sliding a window of the change metric over
   the indicator timeseries.

If `t` (the time vector of `x`), is not provided, it is assumed `t = eachindex(x)`.

Return the output as [`WindowedResults`](@ref) which can be given to
[`significant_transitions`](@ref) to deduce which possible transitions are statistically
significant using a variety of significance tests.
"""
function estimate_transitions(x, config::TransitionsSurrogatesConfig)
    t = eachindex(x)
    return estimate_transitions(t, x, config)
end

function estimate_transitions(config::TransitionsSurrogatesConfig, x, t = eachindex(x))
    # initialize time vectors
    t_indicator = windowmap(config.whichtime, t; width = config.width_ind, stride = config.stride_ind)
    t_change = windowmap(config.whichtime, t_indicator; width = config.width_cha, stride = config.stride_cha)
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
        chametric = one2one ? change_metrics[1] : change_metrics[i]
        z = view(x_indicator, :, i)
        windowmap!(indicator, z, x;
            width = config.width_ind, stride = config.stride_ind
        )
        windowmap!(chametric, view(x_change, :, i), z;
            width = config.width_cha, stride = config.stride_cha
        )
    end

    # put everything together in the output type
    return WindowedResults(
        t, x, t_indicator, x_indicator, t_change, x_change, config
    )
end


"""
    WindowedResults

A struct containing the output of [`estimate_transitions`](@ref) used with
[`WindowedIndicatorMetrics`](@ref).
It can be used for further analysis, visualization,
or given to [`significant_transitions`](@ref).

It has the following fields that the user may access

- `x`: the input timeseries.
- `t`: the time vector of the input timeseries.

- `indicators::Vector{Function}`: indicators used in the processing.
- `x_indicator`, the indicator timeseries (matrix with each column one indicator).
- `t_indicator`, the time vector of the indicator timeseries.

- `change_metrics::Vector{Function}`: change metrics used in the processing.
- `x_change`, the change metric timeseries (matrix with each column one change metric).
- `t_change`, the time vector of the change metric timeseries.
"""
struct TransitionsResults{TT, T<:Real, X<:Real, XX<:AbstractVector{X}, W}
    t::TT # original time vector; most often it is `Base.OneTo`.
    x::XX
    t_indicator::Vector{T}
    x_indicator::Matrix{X}
    t_change::Vector{T}
    x_change::Matrix{X}
    wim::W
end

function Base.show(io::IO, ::MIME"text/plain", res::WindowedResults)
    println(io, "WindowedResults")
    descriptors = [
        "input timeseries" => summary(res.x),
        "indicators" => [nameof(i) for i in res.wim.indicators],
        "indicator (window, stride)" => (res.wim.width_ind, res.wim.stride_ind),
        "change metrics" => [nameof(c) for c in res.wim.change_metrics],
        "change metric (window, stride)" => (res.wim.width_cha, res.wim.stride_cha),
    ]
    padlen = maximum(length(d[1]) for d in descriptors) + 2
    for (desc, val) in descriptors
        println(io, rpad(" $(desc): ", padlen), val)
    end
end