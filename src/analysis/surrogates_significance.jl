"""
    TransitionsSignificance

Supertype used to test for significance in [`significant_transitions`](@ref).
Valid subtypes are:

- [`SurrogatesSignificance`](@ref).
"""
abstract type TransitionsSignificance end


"""
    SurrogatesSignificance <: TransitionsSignificance
    SurrogatesSignificance(; surrogate = RandomFourier(), n = 10_000, tail = :both, rng)

A configuration struct containing instructions on how to test for significance in
the function [`significant_transitions`](@ref) when combined with the output
of [`estimate_transitions`](@ref).

## Description

When used with [`WindowedIndicatorResults`](@ref), significance is estimated as follows:
`n` surrogates from the input timeseries are generated using `surrogate`, which is
any `Surrogate` subtype provided by
[TimeseriesSurrogates.jl](https://juliadynamics.github.io/TimeseriesSurrogates.jl/dev/api/).
For each surrogate, the indicator and then change metric timeseries is extracted.
The values of the surrogate change metrics form a distribution of values (one at each
time point). The value of the original change metric is compared to that of the surrogate
distribution and a p-value is extracted according to the specified `tail`.

the p-value is simply the proportion of surrogate values
that exceed (for `tail = :right`) or subseed (`tail = :left`) the discriminatory
statistic computed from the input data. Use `tail = :left` if
the surrogate data are expected to have higher
discriminatory statistic values. This is the case for statistics that quantify entropy.
For statistics that quantify autocorrelation, use `tail = :right` instead.
For anything else, use the default `tail = :both`.

An iterable of `tail` values can also be given, in which case a specific `tail`
is used for each change metric in [`WindowedIndicatorResults`](@ref).

Keyword `rng = Random.default_rng()` may specify a random number generator for the
surrogates.
"""
Base.@kwdef struct SurrogatesSignificance{S<:Surrogate, T, R} <: TransitionsSignificance
    surrogate::S = RandomFourier()
    n::Int = 10_000
    tail::T = :both
    rng::R = Random.default_rng()
end


"""
    significant_transitions(res::WindowedIndicatorResults, signif::SurrogatesSignificance)

Estimate significant transtions in `res` using the method described by
[`SurrogatesSignificance`](@ref).
Return `pvalues`, a matrix with identical size as `res.x_change`.
It contains the associated p-value timeseries for each change metric (each column).
"""
function significant_transitions(res::WindowedIndicatorResults,
    signif::SurrogatesSignificance)
    (; indicators, change_metrics) = res.wim
    pvalues = similar(res.x_change)
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
        tail = signif.tail isa Symbol ? signif.tail : signif.tail[i]
        c = view(res.x_change, :, i) # change metric timeseries
        # p values for current i
        pval = view(pvalues, :, i)
        indicator_metric_surrogates_loop!(
            indicator, chametric, c, pval, signif.n, sgens,
            indicator_dummys, change_dummys,
            res.wim.width_ind, res.wim.stride_ind, res.wim.width_cha, res.wim.stride_cha, tail
        )
    end
    pvalues ./= signif.n
    return pvalues
end

function indicator_metric_surrogates_loop!(
        indicator, chametric, c, pval, n_surrogates, sgens,
        indicator_dummys, change_dummys,
        width_ind, stride_ind, width_cha, stride_cha, tail
    )

    if tail == :both
        pval_right = zeros(length(pval))
        pval_left = copy(pval_right)
    end

    # parallelized surrogate loop
    Threads.@threads for _ in 1:n_surrogates
        id = Threads.threadid()
        s = sgens[id]()
        change_dummy = change_dummys[id]
        windowmap!(indicator, indicator_dummys[id], s;
            width = width_ind, stride = stride_ind)
        windowmap!(chametric, change_dummy, indicator_dummys[id];
            width = width_cha, stride = stride_cha
        )
        # accumulate for p-value
        if tail == :right
            pval .+= c .< change_dummys[id]
        elseif tail == :left
            pval .+= c .> change_dummy
        elseif tail == :both
            pval_right .+= c .< change_dummy
            pval_left .+= c .> change_dummy
        end
    end
    if tail == :both
        pval .= 2min.(pval_right, pval_left)
    end
end