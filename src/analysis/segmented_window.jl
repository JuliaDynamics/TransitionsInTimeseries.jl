"""
    SegmentedWindowConfig <: IndicatorsChangeConfig
    SegmentedWindowConfig(indicators, change_metrics, tseg_start, tseg_end; kwargs...)

A configuration that can be given to [`estimate_changes`](@ref).
It estimates transitions by estimating indicators and changes in user-defined
window segments as follows:

1. For each segment specified, estimate the corresponding indicator timeseries
   by sliding a window over the input timeseries (within the window segment).
2. For each segment of the indicator timeseries, estimate a scalar change metric by applying
   the change metric over the full segment of the indicator timeseries.d

`tseg_start, tseg_end` are the starts and ends of the window segments
(the window segments may overlap, that's okay).
`indicators, change_metrics` are identical as in [`SlidingWindowConfig`](@ref).

The results output corresponding to `SlidingWindowConfig` is [`SegmentedWindowResults`](@ref).

## Keyword arguments
- `width_ind::Int=100, stride_ind::Int=1, whichtime = midpoint, T = Float64`: keywords
  identical as in [`SlidingWindowConfig`](@ref).
- `min_width_cha::Int=typemax(Int)`: minimal width required to perform the change metric estimation.
  If a segment is not sufficiently long, the change metric is `NaN`.
"""
struct SegmentedWindowConfig{F, G, W<:Function} <: ChangesConfig
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
        min_width_cha = typemax(Int),
        whichtime =  midpoint,
        T = Float64,
    )
    if length(tseg_start) ≠ length(tseg_end)
        throw(ArgumentError("The vectors containing the start and end time of the"*
            " segments must be of equal length."))
    end
    if length(tseg_start) < 1
        throw(ArgumentError("SegmentIndexing requires at least one segment."))
    end
    if count(tseg_start .>= tseg_end) > 1
        throw(ArgumentError("End time of segment must be strictly larger than start time."))
    end
    indicators, change_metrics = sanitycheck_metrics(indicators, change_metrics)
    # Last step: precomputable functions, if any
    indicators = map(f -> precompute(f, 1:T(width_ind)), indicators)

    return SegmentedWindowConfig(
        indicators, change_metrics, tseg_start, tseg_end,
        width_ind, stride_ind, min_width_cha, whichtime,
    )
end

function estimate_changes(config::SegmentedWindowConfig, x, t)
    X, T = eltype(x), eltype(t)
    (; indicators, change_metrics, tseg_start, tseg_end) = config
    n_ind = length(indicators)
    i1, i2 = segment_time(t, tseg_start, tseg_end)

    t_indicator = [T[] for _ in eachindex(tseg_start)]
    x_indicator = [X[;;] for _ in eachindex(tseg_start)]
    t_change = T.(config.tseg_end)
    x_change = fill(X(NaN), length(tseg_start), n_ind)
    one2one = length(change_metrics) == length(indicators)

    dummy_chametric = map(f -> precompute(f, X.(1:2)), change_metrics)
    precomp_change_metrics = [dummy_chametric for _ in eachindex(i1)]

    for k in eachindex(tseg_start)
        tseg = view(t, i1[k]:i2[k])
        xseg = view(x, i1[k]:i2[k])
        t_indicator[k] = windowmap(config.whichtime, tseg; width = config.width_ind,
            stride = config.stride_ind)
        len_ind = length(t_indicator[k])

        # Init with NaN instead of 0 to easily recognise when the segment was too short
        # for the computation to be performed.
        x_indicator[k] = fill(X(NaN), len_ind, n_ind)

        # only analyze if segment long enough to compute metrics
        if len_ind > config.min_width_cha
            precomp_change_metrics[k] = map(f -> precompute(f, t_indicator[k]), change_metrics)
            # Loop over indicators
            for i in 1:n_ind
                indicator = indicators[i]
                chametric = one2one ? precomp_change_metrics[k][i] :
                    precomp_change_metrics[k][1]
                z = view(x_indicator[k], :, i)
                if isnothing(indicator)
                    z = xseg[config.width_ind:config.stride_ind:end]
                else
                    windowmap!(indicator, z, xseg;
                        width = config.width_ind, stride = config.stride_ind
                    )
                end
                if !isnothing(chametric)
                    x_change[k, i] = chametric(z)
                end
            end
        end
    end
    # put everything together in the output type
    return SegmentedWindowResults(t, x, t_indicator, x_indicator, t_change, x_change,
        config, i1, i2, precomp_change_metrics)
end

function segment_time(t::AbstractVector, t1, t2)
    i1 = [searchsortedfirst(t, tt) for tt in t1]
    i2 = [searchsortedlast(t, tt) for tt in t2]
    return i1, i2
end

"""
    SegmentedWindowResults <: ChangesResults

A struct containing the output of [`estimate_changes`](@ref) used with
[`SegmentedWindowConfig`](@ref). It can be used for further analysis, visualization,
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

- `config::SegmentedWindowConfig`, what was used for the analysis.
- `i1::Vector{Int}` indices corresponding to start time of each segment.
- `i2::Vector{Int}` indices corresponding to end time of each segment.
- `precomp_change_metrics` vector containing the precomputed change metrics of each segment.
"""
struct SegmentedWindowResults{TT, T<:Real, X<:Real, XX<:AbstractVector{X},
    W, Z} <: ChangesResults
    t::TT # original time vector; most often it is `Base.OneTo`.
    x::XX
    t_indicator::Vector{Vector{T}}
    x_indicator::Vector{Matrix{X}}
    t_change::Vector{T}
    x_change::Matrix{X}
    config::W
    i1::Vector{Int}
    i2::Vector{Int}
    precomp_change_metrics::Z
end

# Segmented and Sliding results share their show method
function Base.show(io::IO, ::MIME"text/plain", res::Union{SegmentedWindowResults, SlidingWindowResults})
    println(io, nameof(typeof(res)))
    descriptors = [
        "input timeseries" => summary(res.x),
        "indicators" => isnothing(res.config.indicators) ? "nothing" : [nameof(i) for i in res.config.indicators],
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

function show_changemetric(res::SegmentedWindowResults)
    return "change metric (window)" => ([length(t) for t in res.t_indicator])
end
