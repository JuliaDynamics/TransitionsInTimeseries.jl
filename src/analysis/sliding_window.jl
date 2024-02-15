"""
    SlidingWindowConfig <: ChangesConfig
    SlidingWindowConfig(indicators, change_metrics; kwargs...)

A configuration that can be given to [`estimate_indicator_changes`](@ref).
It estimates transitions by a sliding window approach:

1. Estimate the timeseries of an indicator by sliding a window over the input timeseries.
2. Estimate changes of an indicator by sliding a window of the change metric over
   the indicator timeseries.

Both indicators and change metrics are **generic Julia functions** that input an
`x::AbstractVector` and output an `s::Real`. Any function may be given and
see [making custom indicators/change metrics](@ref own_indicator) in the documentation
for more information on possible optimizations.

`indicators` can be a single function or a tuple of indicators.
Similarly, `change_metrics` can be a tuple or a single function.
If tuples, the length of `indicators` and `change_metrics` must match.
This way the analysis
can be efficiently repeated for many indicators and/or change metrics.

The results output corresponding to `SlidingWindowConfig` is [`SlidingWindowResults`](@ref).

Step 1. is skipped if `nothing` is provided as `indicators`, in which case
the change metrics are estimated directly from input data.

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
struct SlidingWindowConfig{F, G, W<:Function} <: ChangesConfig
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
    if !isnothing(indicators)
        indicators = map(f -> precompute(f, 1:T(width_ind)), indicators)
    end
    change_metrics = map(f -> precompute(f, 1:T(width_cha)), change_metrics)

    return SlidingWindowConfig(
        indicators, change_metrics,
        width_ind, stride_ind, width_cha, stride_cha, whichtime,
    )
end

function sanitycheck_metrics(indicators, change_metrics)
    if !(change_metrics isa Tuple)
        change_metrics = (change_metrics,)
    end
    if indicators isa Function
        indicators = (indicators, )
    end
    if !isnothing(indicators) && (length(change_metrics) ≠ length(indicators))
        throw(ArgumentError("The amount of change metrics and indicators must match."))
    end
    return indicators, change_metrics
end

# TODO: This function needs to be split into two as per good practices
# for performant code. One function does the initialization,
# and the other overwrites in place. This makes it easier
# to apply for surrogates as well!

function estimate_indicator_changes(config::SlidingWindowConfig, x, t = eachindex(x))
    (; indicators, change_metrics) = config
    # initialize time vectors
    if isnothing(indicators)
        # Skip indicators if they are nothing
        t_indicator = t
    else
        t_indicator = windowmap(config.whichtime, t; width = config.width_ind,
            stride = config.stride_ind)
    end
    t_change = windowmap(config.whichtime, t_indicator;
        width = config.width_cha, stride = config.stride_cha
    )
    len_ind = length(t_indicator)
    len_change = length(t_change)

    # initialize array containers
    X = eltype(x)
    n_metrics = length(change_metrics)
    x_indicator = zeros(X, len_ind, n_metrics)
    x_change = zeros(X, len_change, n_metrics)

    for i in 1:n_metrics
        # estimate indicator timeseries
        if !isnothing(indicators)
            z = view(x_indicator, :, i)
            windowmap!(indicators[i], z, x;
                width = config.width_ind, stride = config.stride_ind
            )
        else
            # Just equate x here, we're skipping the indicator estimation
            z = x
        end
        # and then apply the change metric
        windowmap!(change_metrics[i], view(x_change, :, i), z;
            width = config.width_cha, stride = config.stride_cha
        )
    end

    # put everything together in the output type
    return SlidingWindowResults(
        t, x, t_indicator, x_indicator, t_change, x_change, config
    )
end

"""
    SlidingWindowResults <: ChangesResults

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

- `config::SlidingWindowConfig`, what was used for the analysis.
"""
struct SlidingWindowResults{TT, T<:Real, X<:Real, XX<:AbstractVector{X}, IT,
    W} <: ChangesResults
    t::TT # original time vector; most often it is `Base.OneTo`.
    x::XX
    t_indicator::IT
    x_indicator::Matrix{X}
    t_change::Vector{T}
    x_change::Matrix{X}
    config::W
end
