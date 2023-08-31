"""
    SignificanceConfig

Supertype used to test for significance in [`estimate_significance`](@ref).
Valid subtypes are:

- [`SurrogatesConfig`](@ref).
"""
abstract type SignificanceConfig end


"""
    SurrogatesConfig <: SignificanceConfig
    SurrogatesConfig(; surrogate = RandomFourier(), n = 1_000, tail = :both, rng)

A configuration struct for significance testing [`estimate_significance`](@ref)
using timeseries surrogates.

## Keyword arguments

- `surromethod = RandomFourier()`: method to generate surrogates
- `n = 1_000`: how many surrogates to generate
- `rng = Random.default_rng()`: random number generator for the surrogates
- `p = 0.05`: threshold for significance of the p-value
- `tail = :both`: tail type used, see below
- `detrend_surro = false`: Boolean specifying whether the surrogates are to be
 detrended by a linear regression. This is sometimes performed in the literature
 but seldomly impacts the results in a significant way, since the surrogates of
 the residual should already be largely free of trend.

## Description

When used with [`WindowedIndicatorResults`](@ref), significance is estimated as follows:
`n` surrogates from the input timeseries are generated using `surromethod`, which is
any `Surrogate` subtype provided by
[TimeseriesSurrogates.jl](https://juliadynamics.github.io/TimeseriesSurrogates.jl/dev/api/).
For each surrogate, the indicator and then change metric timeseries is extracted.
The values of the surrogate change metrics form a distribution of values (one at each time point).
The value of the original change metric is compared to that of the surrogate distribution
and a p-value is extracted according to the specified `tail`.
The p-value is compared with `p` to claim significance, which results in flags.
These results are then stored in [`SurrogatesSignificance`](@ref).

The p-value is simply the proportion of surrogate change metric values
that exceed (for `tail = :right`) or subseed (`tail = :left`) the original change metric
at each given time point. Use `tail = :left` if the surrogate data are expected to have 
higher change metric, discriminatory statistic values. This is the case for statistics 
that quantify entropy. For statistics that quantify autocorrelation, use `tail = :right`
instead. For anything else, use the default `tail = :both`.
An iterable of `tail` values can also be given, in which case a specific `tail`
is used for each change metric in [`WindowedIndicatorResults`](@ref).
"""
struct SurrogatesConfig{S<:Surrogate, T, R} <: SignificanceConfig
    surrogate::S
    n::Int
    tail::T
    rng::R
    p::Float64
    detrend_surro::Bool
end

function SurrogatesConfig(;
    surromethod = RandomFourier(),
    n::Int = 10_000,
    tail = :both,
    rng = Random.default_rng(),
    p = 0.05,
    detrend_surro = false)
    return SurrogatesConfig(surromethod, n, tail, rng, p, detrend_surro)
end

"""
    TransitionsSignificance

Supertype used to store the results of [`estimate_significance`](@ref) based on a
`SignificanceConfig`. Valid subtypes are:

- [`SurrogatesSignificance`](@ref).
"""
abstract type TransitionsSignificance end

"""
    SurrogatesSignificance

A struct containing the results of [`estimate_significance`](@ref) based on a
`SurrogatesConfig`, i.e. the p-values and their associated flags. The former
gives a continuous value between 0 and 1, with 0 denoting an extremely unlikely
and therefore probably significant behavior of the original time series compared
to the ensemble of surrogates. The latter is a Boolean matrix which contains `true`
wherever the p-values are below the threshold defined in [`SurrogatesConfig`](@ref).
Both `signif.p_values` and `signif.flags` have the same size as `res.x_change`, with
`signif::SurrogatesSignificance` and `res::WindowedIndicatorResults`.
"""
mutable struct SurrogatesSignificance <: TransitionsSignificance
    config::SurrogatesConfig
    pvalues::Matrix{Float64}
    flags::BitMatrix
end

"""
    estimate_significance(sigconfig::SignificanceConfig, res::WindowedIndicatorResults)

Estimate the significance of the transition metrics stored in `res` by using the method
described by `sigconfig` and return the results in [`TransitionsSignificance`](@ref).
"""
function estimate_significance(sigconfig::SurrogatesConfig, res::WindowedIndicatorResults)
    (; indicators, change_metrics) = res.wim
    tail = sigconfig.tail
    if !(tail isa Symbol) && length(tail) ≠ length(indicators)
        throw(ArgumentError("Given `tail` must be a symbol or an iterable of same length "*
        "as the input indicators. Got length $(length(tail)) instead of $(length(indicators))."
        ))
    end
    pvalues = zeros(eltype(res.x_change), size(res.x_change)...)

    # Multi-threaded surrogate realization
    seeds = rand(sigconfig.rng, 1:typemax(Int), Threads.nthreads())
    sgens = [surrogenerator(res.x, sigconfig.surrogate, Random.Xoshiro(seed)) for seed in seeds]
    # Dummy vals for surrogate parallelization
    indicator_dummys = [similar(res.x_indicator[:, 1]) for _ in 1:Threads.nthreads()]
    change_dummys = [similar(res.x_change[:, 1]) for _ in 1:Threads.nthreads()]

    for i in 1:size(pvalues, 2) # loop over change metrics
        indicator = indicators[i]
        i_metric = length(change_metrics) == length(indicators) ? i : 1
        chametric = change_metrics[i_metric]
        tai = tail isa Symbol ? tail : tail[i]
        c = view(res.x_change, :, i) # change metric timeseries
        # p values for current i
        pval = view(pvalues, :, i)
        indicator_metric_surrogates_loop!(
            indicator, chametric, c, pval, sigconfig.n, sgens, res.t,
            indicator_dummys, change_dummys, sigconfig.detrend_surro,
            res.wim.width_ind, res.wim.stride_ind, res.wim.width_cha, res.wim.stride_cha, tai
        )
    end
    pvalues ./= sigconfig.n
    return SurrogatesSignificance(sigconfig, pvalues, pvalues .< sigconfig.p)
end

function indicator_metric_surrogates_loop!(
        indicator, chametric, c, pval, n_surrogates, sgens, t,
        indicator_dummys, change_dummys, detrend_surro,
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

        if detrend_surro
            m, p = ridgematrix(t, 0.0) * s
            s .-= (m .* t  .+ p)
        end

        change_dummy = change_dummys[id]
        windowmap!(indicator, indicator_dummys[id], s;
            width = width_ind, stride = stride_ind)
        windowmap!(chametric, change_dummy, indicator_dummys[id];
            width = width_cha, stride = stride_cha)
        # accumulate for p-value
        if tail == :right
            pval .+= c .< change_dummy
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