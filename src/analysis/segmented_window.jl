"""
    SegmentedWindowConfig <: IndicatorsChangeConfig
    SegmentedWindowConfig(indicators, change_metrics, tseg_start, tseg_end; kwargs...)

A configuration that can be given to [`estimate_indicator_changes`](@ref).
It estimates transitions by estimating indicators and changes in user-defined
window segments as follows:

1. For each segment specified, estimate the corresponding indicator timeseries
   by sliding a window over the input timeseries (within the window segment).
2. For each segment of the indicator timeseries, estimate a scalar change metric by applying
   the change metric over the full segment of the indicator timeseries.d

`tseg_start, tseg_end` are the starts and ends of the window segments
(the window segments may overlap, that's okay).
`indicators, change_metrics` are identical as in [`SlidingWindowConfig`](@ref).

## Keyword arguments
-
- `width_ind::Int=100, stride_ind::Int=1, whichtime = midpoint, T = Float64`: keywords
  identical as in [`SlidingWindowConfig`](@ref).
- `min_width_cha::Int=50`: minimal width required to perform the change metric estimation.
  If segment not sufficiently long, return `NaN`.
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

function estimate_indicator_changes(config::SegmentedWindowConfig, x, t)
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

        # Init with NaN instead of 0 to easily recognise when the segment was too short
        # for the computation to be performed.
        x_indicator[k] = fill(NaN, len_ind, n_ind)

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
