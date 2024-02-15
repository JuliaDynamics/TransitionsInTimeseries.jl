# API

## Main analysis functions

```@docs
ChangesConfig
SlidingWindowConfig
SegmentedWindowConfig
estimate_indicator_changes
ChangesResults
SlidingWindowResults
SegmentedWindowResults
```

## Significance testing

```@docs
significant_transitions
TransitionsSignificance
SurrogatesSignificance
ThresholdSignificance
SigmaSignificance
QuantileSignificance
```

## [Indicators](@id indicators)

### Value distribution

```@docs
StatsBase.mean
StatsBase.skewness
StatsBase.kurtosis
```

### Critical Slowing Down

```@docs
StatsBase.var
ar1_whitenoise
```

### Spectrum

```@docs
LowfreqPowerSpectrum
PrecomputedLowfreqPowerSpectrum
```

### Nonlinear dynamics

Indicators that come from nonlinear timeseries analysis typically quantify some entropy-like or complexity measure from the timeseries.
Thousands (_literally_) such measures are provided out of the box by the [ComplexityMeasures.jl](https://juliadynamics.github.io/DynamicalSystemsDocs.jl/complexitymeasures/stable/) package. Given that any of these may be used as an indicator or change metric, we made the decision to not copy-paste any measure here, as it is easy for the user to use any of them.

For example, using the permutation entropy as an indicator is as simple as doing
```julia
using ComplexityMeasures
est = OrdinalPatterns(; m = 3) # order 3
# create a function that given timeseries returns permutation entropy
indicator = x -> entropy_normalized(est, x)
```
and giving the created `indicator` to e.g., [`SlidingWindowConfig`](@ref).

## [Change metrics](@id change_metrics)

### Slope

```@docs
kendalltau
spearman
RidgeRegressionSlope
PrecomputedRidgeRegressionSlope
```

### Value distribution differences

```@docs
difference_of_means
```

## [Make your own indicator/metric!](@id own_indicator)

The only difference between what is an "indicator" and what is a "change metric" is purely conceptual. As far as the code base of TransitionsInTimeseries.jl is concerned, they are both functions `f: x::AbstractVector{Real} -> f(x)::Real`. As a user you may give any such function for an indicator or change metric.

There are situations where you may optimize such a function based on knowledge of input `x` type and length, in which case you may use `PrecomputableFunction`:

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

## Visualization

```@docs
plot_indicator_changes
plot_significance!
plot_changes_significance
```

## Utils

```docs
default_window_width
```