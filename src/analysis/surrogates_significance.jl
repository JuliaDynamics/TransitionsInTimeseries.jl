"""
    TransitionsSignificance

Supertype used to test for significance in [`significant_transitions`](@ref).
Valid subtypes are:

- [`SurrogatesSignificance`](@ref).
"""
abstract type TransitionsSignificance end


"""
    SurrogatesSignificance <: TransitionsSignificance
    SurrogatesSignificance(; surrogate = RandomFourier(), n::Int = 10_000, tail = :both)

A configuration struct containing instructions on how to test for significance in
the function [`significant_transitions`](@ref) when combined with the output
of [`estimate_transitions`](@ref).

## Description

When used with [`WindowedIndicatorResults`](@ref), significance is estimated as follows:
`n` surrogates from the input timeseries are generated using `surrogate`, which is
any `Surrogate` subtype provided by
[TimeseriesSurrogates.jl](https://juliadynamics.github.io/TimeseriesSurrogates.jl/dev/api/).
For each surrogate, the indicator and then change metric is extracted.
The values of the surrogate change metrics form a distribution of values (one at each time point).
The value of the original change metric is compared to that of the surrogate distribution
and a p-value is extracted according to the specified `tail`.

the p-value is simply the proportion of surrogate values
that exceed (for `tail = :right`) or subseed (`tail = :left`) the discriminatory
statistic computed from the input data. Use `tail = :right` if
the surrogate data are expected to have higher
discriminatory statistic values. This is the case for statistics that quantify entropy.
For statistics that quantify autocorrelation, use `tail = :right` instead.
For anything else, use the default `tail = :both`.

An iterable of `tail` values can also be given, in which case a specific `tail`
is used for each change metric in [`WindowedIndicatorResults`](@ref).
"""
Base.@kwdef struct SurrogatesSignificance{S<:Surrogate, T} <: TransitionsSignificance
    surrogate::S = RandomFourier()
    n::Int = 10_000
    tail::T = :both
end
