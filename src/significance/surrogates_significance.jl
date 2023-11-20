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
    # Unpack structs + sanity check
    (; x, x_indicator, x_change) = res
    (; indicators, change_metrics, width_ind, stride_ind, width_cha, stride_cha) = res.config
    (; surrogate, n, tail, rng, p, pvalues) = signif
    n_ind = length(indicators)
    sanitycheck_tail(tail, n_ind)
    
    # Init pvalues
    X = eltype(x_change)
    pvalues = signif.pvalues = zeros(X, size(x_change))
    pvals_right = copy(pvalues)
    pvals_left = copy(pvalues)

    # Multi-threaded surrogate realization
    seeds = rand(rng, 1:typemax(Int), Threads.nthreads())
    sgens = [surrogenerator(x, surrogate, Xoshiro(seed)) for seed in seeds]
    # Dummy vals for surrogate parallelization
    indicator_dummys = [x_indicator[:, 1] for _ in 1:Threads.nthreads()]
    change_dummys = [x_change[:, 1] for _ in 1:Threads.nthreads()]

    Threads.@threads for _ in 1:n

        id = Threads.threadid()
        s = sgens[id]()
        change_dummy = change_dummys[id]

        for i in 1:n_ind
            indicator, change_metric, tai = choose_metrics(indicators, change_metrics,
                tail, i)
            c = view(x_change, :, i)    # change metric timeseries
            pval_right = view(pvals_right, :, i)
            pval_left = view(pvals_left, :, i)
            windowmap!(indicator, indicator_dummys[id], s;
                width = width_ind, stride = stride_ind)
            windowmap!(change_metric, change_dummys[id], indicator_dummys[id];
                width = width_cha, stride = stride_cha)
            accumulate_pvals!(pval_right, pval_left, tai, c, change_dummy)
        end
    end
    choose_pval!(pvalues, pvals_right, pvals_left, tail, n_ind)
    pvalues ./= n
    return pvalues .< p
end

function significant_transitions(res::SegmentedWindowResults, signif::SurrogatesSignificance)
    # Unpack structs + sanity check
    (; x, t_indicator, x_change, i1, i2, precomp_change_metrics) = res
    (; indicators, change_metrics, width_ind, stride_ind, min_width_cha) = res.config
    (; surrogate, n, tail, rng, p) = signif
    n_ind = length(indicators)
    sanitycheck_tail(tail, n_ind)

    # Init pvalues
    X = eltype(x_change)
    pvalues = signif.pvalues = zeros(X, size(x_change))
    pvals_right = copy(pvalues)
    pvals_left = copy(pvalues)

    # Multi-threaded surrogate realization
    seeds = rand(rng, 1:typemax(Int), Threads.nthreads())
    change_dummys = zeros(X, Threads.nthreads())
    xind_length = map(x -> length(x), t_indicator)
    indicator_dummys = Vector{X}[zeros(X, maximum(xind_length)) for _ in
        1:Threads.nthreads()]

    # Loop over segments
    for k in eachindex(i1)

        # If segment too short, return NaN p-value
        if xind_length[k] < min_width_cha
            view(pvalues, k, :) .= X(NaN)
        else
            # Generate surrogates of segment size
            sgens = [surrogenerator(x[i1[k]:i2[k]], surrogate,
                Random.Xoshiro(seed)) for seed in seeds]

            Threads.@threads for _ in 1:n
                id = Threads.threadid()
                s = sgens[id]()
                for i in eachindex(indicators)
                    indicator, chametric, tai = choose_metrics(indicators,
                        precomp_change_metrics[k], tail, i)
                    c = x_change[k, i]
                    indicator_dummy = view(indicator_dummys[id], 1:xind_length[k])
                    windowmap!(indicator, indicator_dummy, s;
                        width = width_ind, stride = stride_ind)
                    change_dummys[id] = chametric(indicator_dummy)
                    accumulate_pvals!(pvals_right, pvals_left, tai, c, change_dummys[id],
                        k, i)
                    choose_pval!(pvalues, pvals_right, pvals_left, tail, k, i)
                end
            end
        end
    end
    pvalues ./= n
    return pvalues .< p
end

function sanitycheck_tail(tail, n_ind)
    if !(tail isa Symbol) && length(tail) ≠ n_ind
        throw(ArgumentError("Given `tail` must be a symbol or an iterable of same length "*
        "as the input indicators. Got length $(length(tail)) instead of $n_ind."
        ))
    end
end

function choose_metrics(indicators::Tuple{Vararg}, change_metrics, tail, i::Int)
    i_metric = length(change_metrics) == length(indicators) ? i : 1
    tai = tail isa Symbol ? tail : tail[i]
    return indicators[i], change_metrics[i_metric], tai
end

function accumulate_pvals!(pval_right, pval_left, tail, c, change_dummy)
    if tail == :both || tail == :right
        pval_right .+= c .< change_dummy
    end
    if tail == :both || tail == :left
        pval_left .+= c .> change_dummy
    end
end

function accumulate_pvals!(pvals_right, pvals_left, tail, c, change_dummy, k, i)
    if tail == :both || tail == :right
        pvals_right[k, i] += c < change_dummy
    end
    if tail == :both || tail == :left
        pvals_left[k, i] += c > change_dummy
    end
end

# loop vectorial choose_pval! over indices
function choose_pval!(pvals, pvals_right, pvals_left, tail, n_ind)
    for i in 1:n_ind
        pval = view(pvals, :, i)
        pval_right = view(pvals_right, :, i)
        pval_left = view(pvals_left, :, i)
        choose_pval!(pval, pval_right, pval_left, tail)
    end
end

# vectorial choose_pval!
function choose_pval!(pval, pval_right, pval_left, tail)
    if tail == :both
        pval .= 2min.(pval_right, pval_left)
    elseif tail == :right
        pval .= pval_right
    elseif tail == :left
        pval .= pval_left
    end
end

# scalar choose_pval!
function choose_pval!(pvals, pvals_right, pvals_left, tail, k, i)
    if tail == :both
        pvals[k, i] = 2min(pvals_right[k, i], pvals_left[k, i])
    elseif tail == :right
        pvals[k, i] = pvals_right[k, i]
    elseif tail == :left
        pvals[k, i] = pvals_left[k, i]
    end
end