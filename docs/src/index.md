# TransitionIndicators.jl

![TransitionIndicators.jl](assets/logo.gif)

```@docs
TransitionIndicators
```

!!! info "Star us on GitHub!"
    If you have found this package useful, please consider starring it on [GitHub](https://github.com/JuliaDynamics/TransitionIndicators.jl).
    This gives us an accurate lower bound of the (satisfied) user count.

## [Content] (@id content)

Multi-stability is an ubiquitous feature of real-world dynamical systems, as for instance in epidemiology, sociology and climate science. The transition between two stability regimes can be abrupt and of large magnitude, although its trigger can be gradual and of small amplitude. The repercusions of such transitions can be hard to mitigate. Finding mathematical tools to predict their occurence solely based on data has therefore been an active field of research over the last decades. Hereby, many terms were employed such as *early warning signals*, *resilience indicators*, *regime-shift identifiers*, *change-point detection* and *transition indicators*. `TransitionIndicators.jl` sticks to the latter terminology and provides an interface that:

- Allows a fast computation of common transition indicators with a couple of lines, as demonstrated in the [example section](@ref example).
- Makes the surrogate analysis to test for significance under the hub.
- Can be easily extended by any user without touching the source code.
- Reduces the programming overhead for any researcher willing to benchmark new methods.
- Eases the reproducibility thanks to a clear syntax and a simple installation.
- Increases trustworthiness thanks to a large test suite.

!!! info "Similar projects"
    An R [toolbox](https://www.early-warning-signals.org/?page_id=42) and a Python [library](https://pypi.org/project/ewstools/) already exist. However, we believe that they are difficult to extend for the user. Furthermore, they do not offer a native performant code, as here allowed by the use of Julia.


## [Under the hub] (@id workflow)

Computing transition indicators is schmetically represented in the plot below and essentially consists of:
1. Detrending the time series to get the fluctuations around the tracked attractor.
1. Estimating the time series of an indicator by sliding a window over the fluctuation time-series.
2. Computing the evolution of the indicator by sliding a window over its time series.
3. Generating many surrogates that preserve important statistical properties of the originial fluctuation time-series.
4. Performing step 2 and 3 for the surrogate time series to check whether the indicator evolution of the original time series shows a significant feature (trend or jump).

![Schematic representation of what is happening under the hub.](assets/workflow.svg)

## [Example] (@id example)

Let us generate data from a bi-stable nonlinear model subject to noise and to a gradual change of the forcing that leads to a transition. Furthermore, we also generate data from a linear model, which is by definition monostable and therefore incapable of transitionning. This is done to control the rate of false positives. Generating (and plotting) such prototypical data to test some indicators can be done by simply running:

```@example MAIN
using TransitionIndicators
using CairoMakie

t, x_linear, x_nlinear = generate_test_data()
fig = Figure()
axs = [Axis(fig[i, j], xticklabelsvisible = i == 5 ? true : false) for i in 1:5, j in 1:2]
[xlims!(axs[i, j], (0, 50)) for i in 1:5, j in 1:2]
lines!(axs[1, 1], t, x_linear)
lines!(axs[1, 2], t, x_nlinear)
fig
```

Finding transition indicators focuses on analyzing the fluctuations of the system around the tracked attractor. Therefore, a detrending step is common and can be performed by the user, for instance with the help of a wide range of `julia` packages. We particularly recommend the use of [`Loess.jl`](https://github.com/JuliaStats/Loess.jl), an implementation of the Locally Estimated Scatterplot Smoothing (Loess), which has become one of the standard chance for analyzing fluctuations around attractors.

```@example MAIN
using Loess

xlin_loess = predict(loess(t, x_linear, span = 0.1), t)
xnlin_loess = predict(loess(t, x_nlinear, span = 0.1), t)
lines!(axs[1, 1], t, xlin_loess)
lines!(axs[1, 2], t, xnlin_loess)

lin_fluctuations = x_linear - xlin_loess
nlin_fluctuations = x_nlinear - xnlin_loess
lines!(axs[2, 1], t, lin_fluctuations)
lines!(axs[2, 2], t, nlin_fluctuations)

fig
```

Performing the computation of some transition indicators can now be easily done by defining some parameters and providing a vector of indicators and of evolution metrics:

```@example MAIN
p = init_metaanalysis_params(
    wv_indicator_width = 60,
    wv_indicator_stride = 2,
    wv_evolution_width = 30,
    wv_evolution_stride = 2,
)
indicators = [variance, ar1_whitenoise]
evolution_metrics = precomputed_ridge_slope(p)

results_lin = analyze_indicators(t, lin_fluctuations, indicators, evolution_metrics, p)
results_nonlin = analyze_indicators(t, nlin_fluctuations, indicators, evolution_metrics, p)

results = [results_lin, results_nonlin]
indicator_number = 1
for j in eachindex(results)
    lines!(axs[3, j], results[j].t_indicator, results[j].X_indicator[1, :, indicator_number])
    lines!(axs[4, j], results[j].t_evolution, results[j].X_evolution[1, :, indicator_number])
end
fig
```

Once we obtained the results, we would like to know whether the increase of variance and AR1-regression coefficient are significant:

```@example MAIN
sig_lin = measure_significances(results_lin, normalized_confidence_intervall)
sig_nonlin = measure_significances(results_nonlin, normalized_confidence_intervall)
significances = [sig_lin, sig_nonlin]
for i in eachindex(indicators), j in eachindex(results)
    lines!(axs[5, j], results[j].t_evolution, significances[j][1, :, i])
end
fig
```

### Copy-pastable code

```julia
using TransitionIndicators
using CairoMakie
using Loess

t, x_linear, x_nlinear = generate_test_data()
fig = Figure()
axs = [Axis(fig[i, j], xticklabelsvisible = i == 5 ? true : false) for i in 1:5, j in 1:2]
[xlims!(axs[i, j], (0, 50)) for i in 1:5, j in 1:2]
lines!(axs[1, 1], t, x_linear)
lines!(axs[1, 2], t, x_nlinear)

xlin_loess = predict(loess(t, x_linear, span = 0.1), t)
xnlin_loess = predict(loess(t, x_nlinear, span = 0.1), t)
lines!(axs[1, 1], t, xlin_loess)
lines!(axs[1, 2], t, xnlin_loess)
lin_fluctuations = x_linear - xlin_loess
nlin_fluctuations = x_nlinear - xnlin_loess
lines!(axs[2, 1], t, lin_fluctuations)
lines!(axs[2, 2], t, nlin_fluctuations)

p = init_metaanalysis_params(
    wv_indicator_width = 60,
    wv_indicator_stride = 2,
    wv_evolution_width = 30,
    wv_evolution_stride = 2,
)
indicators = [variance, ar1_whitenoise]
evolution_metrics = precomputed_ridge_slope(p)
results_lin = analyze_indicators(t, lin_fluctuations, indicators, evolution_metrics, p)
results_nonlin = analyze_indicators(t, nlin_fluctuations, indicators, evolution_metrics, p)
results = [results_lin, results_nonlin]
indicator_number = 1
for j in eachindex(results)
    lines!(axs[3, j], results[j].t_indicator, results[j].X_indicator[1, :, indicator_number])
    lines!(axs[4, j], results[j].t_evolution, results[j].X_evolution[1, :, indicator_number])
end

sig_lin = measure_significances(results_lin, normalized_confidence_intervall)
sig_nonlin = measure_significances(results_nonlin, normalized_confidence_intervall)
significances = [sig_lin, sig_nonlin]
for i in eachindex(indicators), j in eachindex(results)
    lines!(axs[5, j], results[j].t_evolution, significances[j][1, :, i])
end
fig
```

## API

### Data generation

```@docs
generate_test_data
```

### High-level interface
```@docs
IndicatorEvolutionResults
init_metaanalysis_params
analyze_indicators
```

### Indicators
```@docs
ar1_whitenoise
variance
krt
skw
```

### Evolution metrics
```@docs
precomputed_ridge_slope
kendalltau
spearman
```

### Surrogates

For the surrogate generation, you can use any technique defined in [TimeseriesSurrogates.jl](https://juliadynamics.github.io/TimeseriesSurrogates.jl/stable/#Surrogate-methods).

### Significance metrics

```@docs
measure_significance
measure_significances
normalized_confidence_intervall
which_percentile
normalized_percentile_distance
```

