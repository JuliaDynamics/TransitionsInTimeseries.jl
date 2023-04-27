# API

## Load data

```@docs
load_linear_vs_doublewell()
```

## Params
```@docs
IndicatorsParams
ChangeMetricsParams
init_metrics
```

## Indicators analysis
```@docs
IndicatorsConfig
SignificanceConfig
indicators_analysis
IndicatorsResults
```

## [Indicators] (@id indicator_functions)
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

For the surrogate generation, you can use any subtype of `Surrogate` defined in [Timeseriessurrogates.jl](https://juliadynamics.github.io/Timeseriessurrogates.jl/stable/#surrogate-methods).

## Sliding windows
```@docs
WindowViewer
windowmap
windowmap!
slidebracket
```