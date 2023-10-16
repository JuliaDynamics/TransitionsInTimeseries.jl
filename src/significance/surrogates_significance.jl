"""
    SurrogatesSignificance <: TransitionsSignificance
    SurrogatesSignificance(; surrogate = RandomFourier(), n = 10_000, tail = :both, rng)

A configuration struct for significance testing [`significant_transitions`](@ref)
using timeseries surrogates.

## Keyword arguments

- `surromethod = RandomFourier()`: method to generate surrogates
- `n = 1000`: how many surrogates to generate
- `rng = Random.default_rng()`: random number generator for the surrogates
- `p = 0.05`: threshold for significance of the p-value
- `tail = :both`: tail type used, see below

## Description

When used with [`IndicatorsChangesResults`](@ref), significance is estimated as follows:
`n` surrogates from the input timeseries are generated using `surromethod`, which is
any `Surrogate` subtype provided by
[TimeseriesSurrogates.jl](https://juliadynamics.github.io/TimeseriesSurrogates.jl/dev/api/).
For each surrogate, the indicator and then change metric timeseries is extracted.
The values of the surrogate change metrics form a distribution of values (one at each time point).
The value of the original change metric is compared to that of the surrogate distribution
and a p-value is extracted according to the specified `tail`.
The p-value is compared with `p` to claim significance.
After using `SurrogatesSignificance`, you may access the full p-values before thresholding
in the field `.pvalues` (to e.g., threshold with different `p`).

The p-value is simply the proportion of surrogate change metric values
that exceed (for `tail = :right`) or subseed (`tail = :left`) the original change metric
at each given time point. Use `tail = :left` if the surrogate data are expected to have
higher change metric, discriminatory statistic values. This is the case for statistics
that quantify entropy. For statistics that quantify autocorrelation, use `tail = :right`
instead. For anything else, use the default `tail = :both`.
An iterable of `tail` values can also be given, in which case a specific `tail`
is used for each change metric in [`IndicatorsChangesResults`](@ref).

Note that the raw p-values can be accessed in the field `.pvalues`, after calling the
[`significant_transitions`](@ref) function with `SurrogatesSignificance`, in case you wish
to obtain a different threshold of the p-values.
"""
mutable struct SurrogatesSignificance{S<:Surrogate, T, R} <: TransitionsSignificance
    surrogate::S
    n::Int
    tail::T
    rng::R
    p::Float64
    pvalues::Matrix{Float64}
end

function SurrogatesSignificance(;
    surromethod = RandomFourier(),
    n::Int = 10_000,
    tail = :both,
    rng = Random.default_rng(), p = 0.05)
    return SurrogatesSignificance(surromethod, n, tail, rng, p, zeros(1,1))
end

function significant_transitions(res::SlidingWindowResults, signif::SurrogatesSignificance)
    (; indicators, change_metrics) = res.config
    tail = signif.tail
    if !(tail isa Symbol) && length(tail) ≠ length(indicators)
        throw(ArgumentError("Given `tail` must be a symbol or an iterable of same length "*
        "as the input indicators. Got length $(length(tail)) instead of $(length(indicators))."
        ))
    end
    pvalues = signif.pvalues = similar(res.x_change)
    # Multi-threaded surrogate realization
    seeds = rand(signif.rng, 1:typemax(Int), Threads.nthreads())
    sgens = [surrogenerator(res.x, signif.surrogate, Random.Xoshiro(seed)) for seed in seeds]
    # Dummy vals for surrogate parallelization
    indicator_dummys = [res.x_indicator[:, 1] for _ in 1:Threads.nthreads()]
    change_dummys = [res.x_change[:, 1] for _ in 1:Threads.nthreads()]
    for i in 1:size(pvalues, 2) # loop over change metrics
        indicator = indicators[i]
        i_metric = length(change_metrics) == length(indicators) ? i : 1
        chametric = change_metrics[i_metric]
        tai = tail isa Symbol ? tail : tail[i]
        c = view(res.x_change, :, i) # change metric timeseries
        # p values for current i
        pval = view(pvalues, :, i)
        sliding_surrogates_loop!(
            indicator, chametric, c, pval, signif.n, sgens,
            indicator_dummys, change_dummys,
            res.config.width_ind, res.config.stride_ind, res.config.width_cha, res.config.stride_cha, tai
        )
    end
    pvalues ./= signif.n
    return pvalues .< signif.p
end

function significant_transitions(res::SegmentedWindowResults, signif::SurrogatesSignificance)
    # Unpack and sanity checks
    X = eltype(res.x_change)
    (; indicators, change_metrics, tseg_start, tseg_end) = res.config
    tail = signif.tail
    if !(tail isa Symbol) && length(tail) ≠ length(indicators)
        throw(ArgumentError("Given `tail` must be a symbol or an iterable of same length "*
        "as the input indicators. Got length $(length(tail)) instead of $(length(indicators))."
        ))
    end
    # Multi-threaded surrogate realization
    seeds = rand(signif.rng, 1:typemax(Int), Threads.nthreads())
    change_dummys = [X(0) for _ in 1:Threads.nthreads()]
    pvalues = signif.pvalues = zeros(X, size(res.x_change)...)

    for k in eachindex(tseg_start)  # loop over segments
        # If segment too short, return inf p-value
        if length(res.x_indicator[k][:, 1]) < res.config.min_width_cha
            pvalues[k, :] .= Inf
        else
            tseg, xseg = segment(res.t, res.x, tseg_start[k], tseg_end[k])
            sgens = [surrogenerator(xseg, signif.surrogate, Random.Xoshiro(seed))
                for seed in seeds]
            # Dummy vals for surrogate parallelization
            indicator_dummys = [res.x_indicator[k][:, 1] for _ in 1:Threads.nthreads()]
            for i in eachindex(indicators)
                indicator = indicators[i]
                i_metric = length(change_metrics) == length(indicators) ? i : 1
                chametric = change_metrics[i_metric]

                # Precomputation: include time vector for metrics that require it!
                if chametric isa PrecomputableFunction
                    chametric = precompute(chametric, res.t_indicator[k])
                end

                tai = tail isa Symbol ? tail : tail[i]
                c = res.x_change[k, i]  # change metric
                pval = view(pvalues, k, i)
                # p values for current i
                segmented_surrogates_loop!(
                    indicator, chametric, c, pval, signif.n, sgens,
                    indicator_dummys, change_dummys,
                    res.config.width_ind, res.config.stride_ind, tai)
            end
        end
    end
    pvalues ./= signif.n
    return pvalues .< signif.p
end

function sliding_surrogates_loop!(
        indicator, chametric, c, pval, n_surrogates, sgens,
        indicator_dummys, change_dummys,
        width_ind, stride_ind, width_cha, stride_cha, tail
    )
    pval_right, pval_left = init_pvals(pval)
    # parallelized surrogate loop
    Threads.@threads for _ in 1:n_surrogates
        id = Threads.threadid()
        s = sgens[id]()
        change_dummy = change_dummys[id]
        windowmap!(indicator, indicator_dummys[id], s;
            width = width_ind, stride = stride_ind)
        windowmap!(chametric, change_dummy, indicator_dummys[id];
            width = width_cha, stride = stride_cha)
        accumulate_pvals!(pval_right, pval_left, tail, c, change_dummy)
    end
    choose_pval!(pval, pval_right, pval_left, tail)
end

function segmented_surrogates_loop!(
    indicator, chametric, c, pval, n_surrogates, sgens,
    indicator_dummys, change_dummys, width_ind, stride_ind, tail
)
    pval_right, pval_left = init_pvals(pval)
    # parallelized surrogate loop
    Threads.@threads for _ in 1:n_surrogates
        id = Threads.threadid()
        s = sgens[id]()
        change_dummy = change_dummys[id]
        windowmap!(indicator, indicator_dummys[id], s;
            width = width_ind, stride = stride_ind)
        change_dummy = chametric(indicator_dummys[id])
        accumulate_pvals!(pval_right, pval_left, tail, c, change_dummy)
    end
    choose_pval!(pval, pval_right, pval_left, tail)
end

function init_pvals(pval)
    pval_right = zeros(length(pval))
    pval_left = copy(pval_right)
    return pval_right, pval_left
end

function accumulate_pvals!(pval_right, pval_left, tail, c, change_dummy)
    if tail == :both || tail == :right
        pval_right .+= c .< change_dummy
    end
    if tail == :both || tail == :left
        pval_left .+= c .> change_dummy
    end
end

function choose_pval!(pval, pval_right, pval_left, tail)
    if tail == :both
        println(pval_right, pval_left)
        pval .= 2min.(pval_right, pval_left)
    elseif tail == :right
        pval .= pval_right
    elseif tail == :left
        pval .= pval_left
    end
end
