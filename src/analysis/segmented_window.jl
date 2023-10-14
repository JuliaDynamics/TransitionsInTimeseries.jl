
"""
    SegmentedWindowConfig(indicators, change_metrics, tseg_start, tseg_end;
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
struct SegmentedWindowConfig{F, G, W<:Function} <: IndicatorsChangesConfig
    indicators::F
    change_metrics::G
    tseg_start::Vector
    tseg_end::Vector
    width_ind::Int
    stride_ind::Int
    min_width_cha::Int
    whichtime::W
end

function SegmentedWindowConfig(
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

    return SegmentedWindowConfig(
        indicators, change_metrics, tseg_start, tseg_end,
        width_ind, stride_ind, min_width_cha, whichtime,
    )
end
