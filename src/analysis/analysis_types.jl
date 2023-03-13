struct IndicatorsInput{X<:AbstractVector, F<:Function, W}
    indicators::Vector{<:F}
    window_width::W
    window_stride::W
end

FF = Union{Function, Vector{<:Function}}

function IndicatorsInput(f::FF; kwargs...) where {X<:AbstractVector}
    if f isa Function
        ff = [f]
    else
        ff = f
    end
    return IndicatorsInput{X, typeof(ff), typeof(kwargs)}(x, ff, kwargs)
end


"""

    SignificanceHyperParams(kwargs...)

Initialize a `SignificanceHyperParams` struct that dictates how a significance
meta-analysis will be done on transition indicators.

## Keyword arguments
- `n_surrogates::Int = 10_000`: how many surrogates to create.
- `surrogate_method::S = RandomFourier()`: what method to use to create the surrogates.
  Any `Surrogate` subtype from [TimeseriesSurrogates.jl](
    https://juliadynamics.github.io/TimeseriesSurrogates.jl/stable/#Surrogate-methods
  ) is valid.
- `rng::AbstractRNG = Random`.default_rng()`: a random number generator for the surrogates.
  See the [Julia manual](https://docs.julialang.org/en/v1/stdlib/Random/#Random-Numbers) for more.
- `wv_indicator_width::Int = 100`,
- `wv_indicator_stride::Int = 5`,
- `wv_evolution_width::Int = 50`,
- `wv_evolution_stride::Int = 5`,
"""
function SignificanceHyperParams(;
    n_surrogates::Int = 10_000,
    surrogate_method::S = RandomFourier(),
    rng::AbstractRNG = Random.default_rng(),
    wv_indicator_width::Int = 100,
    wv_indicator_stride::Int = 5,
    wv_evolution_width::Int = 50,
    wv_evolution_stride::Int = 5,
) where {S<:Surrogate}
    return SignificanceHyperParams(
        n_surrogates, surrogate_method, rng,
        wv_indicator_width, wv_indicator_stride,
        wv_evolution_width, wv_evolution_stride,
    )
end



function analyze_indicators(x, indicators, significance_params)



end



"""
    IndicatorsResults

A struct containing the input and output of the main computational part of
TransitionIndicators.jl. It can be ginen to the [`indicator_significance`](@ref) function.

It has the following fields:

## Input

- `t`: the time vector of the input timeseries
- `x`: the input timeseries
- `indicators::Vector{Function}`: indicators used in the processing


- `t_indicator`, the time vector of `X_indicator`. dims = `nt_indicator`.
- `X_indicator`, the indicator time series. dims = 1 x `nt_indicator` x `ni`.
- `t_evolution`, the time vector of `X_evolution` and `S_evolution`. dims = `nt_evolution`.
- `X_evolution`, the time series of the evolution metric. dims = 1 x `nt_evolution` x `ni`.
- `S_evolution`, the time series of the evolution metric computed on the surrogates. dims = `ns` x `nt_evolution` x `ni`.

with:
- `ni`: number of indicators to compute
- `nt_indicator`: number of data points in the indicator time series
- `nt_evolution`: number of data points in the evolution metric time series
- `ns`: number of surrogates
"""
struct IndicatorEvolutionResults{T<:AbstractFloat}
    t_indicator::Vector{T}
    X_indicator::Array{T, 3}
    t_evolution::Vector{T}
    X_evolution::Array{T, 3}
    S_evolution::Array{T, 3}
end

struct SignificanceHyperParams{S<:Surrogate, R<:AbstractRNG}
    n_surrogates::Int
    surrogate_method::S
    rng::R
    wv_indicator_width::Int
    wv_indicator_stride::Int
    wv_evolution_width::Int
    wv_evolution_stride::Int
end

# TODO: init struct based on dimensions of input time-series

#####################################################
# Change-point metrics
#####################################################