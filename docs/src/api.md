# API

## Load data

```@docs
load_linear_vs_doublewell()
```

## High-level interface
```@docs
IndicatorsConfig
SignificanceConfig
indicators_analysis
indicators_significance
```

## [Indicators] (@id indicator_functions)
```@docs
ar1_whitenoise
Statistics.mean(::Any)
Statistics.var(::AbstractArray)
StatsBase.skewness
StatsBase.kurtosis
```

## Change metrics
```@docs
RidgeRegression
kendalltau
spearman
```

## Surrogates

For the surrogate generation, you can use any subtype of `Surrogate` defined in [Timeseriessurrogates.jl](https://juliadynamics.github.io/Timeseriessurrogates.jl/stable/#surrogate-methods).


## Significance metrics

```@docs
threshold_indicators
measure_significance
confidence_interval
normalized_percentile
```

## Sliding windows
```@docs
WindowViewer
windowmap
windowmap!
```