# API

## Load data

```@docs
load_linear_vs_doublewell()
```

## Indicators analysis
```@docs
precompute_metrics
IndicatorsConfig
SignificanceConfig
indicators_analysis
IndicatorsResults
```

## [Indicators](@id Indicators)
```@docs
Statistics.mean(::Any)
Statistics.var(::AbstractArray)
StatsBase.skewness
StatsBase.kurtosis
ar1_whitenoise
LowfreqPowerSpectrum
PermutationEntropy
```

## Change metrics
```@docs
kendalltau
spearman
RidgeRegressionSlope
```

## Surrogates

For the surrogate generation, you can use any subtype of `Surrogate` defined in [Timeseriessurrogates.jl](https://juliadynamics.github.io/TimeseriesSurrogates.jl/v1.0/#Surrogate-methods-1).

## Sliding windows
```@docs
WindowViewer
windowmap
windowmap!
```