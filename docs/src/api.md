# API

## Main analysis functions

```@docs
WindowedIndicatorConfig
transition_metrics
WindowedIndicatorResults
```

## Significance testing

```@docs
estimate_significance
SignificanceConfig
SurrogatesConfig
TransitionsSignificance
SurrogatesSignificance
QuantileSignificance
```

## [Indicators](@id indicators)

### Value distribution

```@docs
Statistics.mean(::Any)
StatsBase.skewness
StatsBase.kurtosis
```

### Critical Slowing Down

```@docs
Statistics.var(::AbstractArray)
ar1_whitenoise
```

### Spectrum

```@docs
LowfreqPowerSpectrum
```

### Nonlinear dynamics

Indicators that come from nonlinear timeseries analysis and quantify some entropy-based dynamic quantity in the timeseries. They are provided by the [ComplexityMeasures.jl](https://juliadynamics.github.io/ComplexityMeasures.jl/stable/) package, that lists 100s of possible such indicators. Here we only provide an indicator out of the box for the permutation entropy, but
building something similar is trivial:
```julia
function permutation_entropy(; m = 3, τ = 1)
    est = SymbolicPermutation(; m, τ)
    return x -> entropy_normalized(est, x)
end
```

```@docs
permutation_entropy
entropy
```

## [Change metrics](@id change_metrics)

### Slope

```@docs
kendalltau
spearman
RidgeRegressionSlope
```

### Value distribution differences

```@docs
difference_of_means
```

## [Make your own indicator/metric!](@id own_indicator)

The only difference between what is an "indicator" and what is a "change metric" is purely conceptual. As far as the code base of TransitionsInTimeseries.jl is concerned, they are both functions `f: x::AbstractVector{Real} -> f(x)::Real`. As a user you may give any such function for an indicator or change metric.

There are situations where you may optimize such a function based on knowledge of input `x` type and length.

TODO: Here explain how to use precomputable functions


```@docs
PrecomputableFunction
precompute
```

## Surrogates

For the surrogate generation, you can use any subtype of `Surrogate` defined in [Timeseriessurrogates.jl](https://juliadynamics.github.io/TimeseriesSurrogates.jl/v1.0/#Surrogate-methods-1).

## Sliding windows
```@docs
WindowViewer
windowmap
windowmap!
```

## Load data

```@docs
load_linear_vs_doublewell()
```
