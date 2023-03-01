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
- Makes the surrogate analysis to test for significance [under the hub](@ref workflow).
- Can be easily extended by any user without touching the source code.
- Reduces the programming overhead for any researcher willing to benchmark new methods.
- Eases the reproducibility thanks to a clear syntax and a simple installation.
- Increases trustworthiness thanks to a large test suite.

!!! info "Similar projects"
    An R [toolbox](https://www.early-warning-signals.org/?page_id=42) and a Python [library](https://pypi.org/project/ewstools/) already exist. However, we believe that they are difficult to extend for the user. Furthermore, they do not offer a native performant code, as here allowed by the use of Julia.

## [Approaches] (@id approaches)

Over the last decades, research on transition indicators has largely focused on Critical Slowing Down (CSD). CSD is observed when a system with continuous right-hand side approaches a bifurcation and consists in a resilience loss of the system. For instance this can be diagnosed by an increase of the variance and the AR1-regression coefficient, as demonstrated in the [example section](@ref example). However, we emphasize that this is one out of many possible approaches for obtaining transition indicators. Recent work has explored new approaches relying on nonlinear dynamics or machine learning. `TransitionIndicators.jl` is designed to allow these cutting-edge methods and foster the development of new ones.

## [Under the hub] (@id workflow)

Computing transition indicators is schmetically represented in the plot below and essentially consists of:
1. Detrending the time series to get the fluctuations around the tracked attractor.
1. Estimating the time series of an indicator by sliding a window over the fluctuation time-series.
2. Computing the evolution of the indicator by sliding a window over its time series.
3. Generating many surrogates that preserve important statistical properties of the originial fluctuation time-series.
4. Performing step 2 and 3 for the surrogate time series to check whether the indicator evolution of the original time series shows a significant feature (trend or jump).

![Schematic representation of what is happening under the hub.](assets/workflow.svg)

!!! info "Adaptability"
    Some approaches mentioned [above](@ref approaches) do not require all of these steps. For instance, a deep-learning classifier might not need the generation of surrogates, as an equivalent step might have been performed at training time. However, note that each of the listed steps can be easily by-passed. `TransitionIndicators.jl` therefore provides a general framework that is equally well suited for specific cases.

## [Example] (@id example)

Let us load data from a bi-stable nonlinear model subject to noise and to a gradual change of the forcing that leads to a transition. Furthermore, we also load data from a linear model, which is by definition monostable and therefore incapable of transitionning. This is done to control the rate of false positives, a common problem that can emerge when looking for tranisiton indicators. Loading such prototypical data to test some indicators can be done by simply running:

```@example MAIN
using TransitionIndicators

t, x_linear, x_nlinear = load_test_data()
X = [x_linear, x_nlinear]
```

Indicating a transition can now be easily done by defining some indicators and evolution metrics to compute based on some parameters initialized with `init_metaanalysis_params`. The detrending step is here performed by simply building the difference to the previous data point. After this, obtaining the time at which the evolution metric of the indicators is significant is a matter of three lines. Furthermore, let us visualize the results by using [`Makie.jl`](https://docs.makie.org/stable/):

```@example MAIN
using CairoMakie

# Initialize the indicator analysis
p = init_metaanalysis_params(
    wv_indicator_width = 400,
    wv_indicator_stride = 1,
    wv_evolution_width = 20,
    wv_evolution_stride = 1,
)
indicators = [variance, ar1_whitenoise]
evolution_metrics = precomputed_ridge_slope(p)

# Initialize figure and axes
fig = Figure()
nrows, ncols = 5, 2
axs = [Axis(
    fig[i, j],
    xticklabelsvisible = i == 5 ? true : false,
) for i in 1:nrows, j in 1:ncols]
[xlims!(axs[i, j], extrema(t)) for i in 1:nrows, j in 1:ncols]

# Choose which indicator (evolution) to plot in row 3 (and 4)
# Here we choose the first one provided in `indicators`
ind_idx = 1

for j in eachindex(X)
    # Indicator analysis
    x = X[j]
    fluctuation = diff(x)
    result = analyze_indicators(t[2:end], fluctuation, indicators, evolution_metrics, p)
    significance = measure_significances(result, normalized_confidence_intervall)
    t_indicator = threshold_indicators(result.t_evolution, significance)
    
    # Plot results
    lines!(axs[1, j], t, x)
    lines!(axs[2, j], t[2:end], fluctuation)
    lines!(axs[3, j], result.t_indicator, result.X_indicator[1, :,ind_idx])
    lines!(axs[4, j], result.t_evolution, result.X_evolution[1, :, ind_idx])
    lines!(axs[5, j], result.t_evolution, significance[1, :, 1])
    lines!(axs[5, j], result.t_evolution, significance[1, :, 2])
    [vlines!(axs[i, j], t_indicator, linestyle = :dash, color = :red) for i in [1, 5]]
end
fig
```

Here we see that a transition is correctly forecasted for the double-fold system, whereas none is predicted for the linear system. Note that you can abitrarily extend the vector specifying which indicators should be computed - either with functions implemented in `TransitionIndicators.jl` and listed [here](@ref indicator_functions), or with user-defined functions that respect the vector-input/scalar-ouput structure assumed. The same holds for the evolution and significance metrics.

!!! info "Detrending in Julia"
    Detrending can be performed in many ways and therefore remains an external step of `TransitionIndicators.jl`. A wide range of `Julia` packages exists to perform smoothing such as [`Loess.jl`](https://github.com/JuliaStats/Loess.jl), [`DSP.jl`](https://docs.juliadsp.org/latest/contents/) or [`RollingFunctions.jl`](https://jeffreysarnoff.github.io/RollingFunctions.jl/dev/). The detrending step then simply consists of subtracting the smoothed signal from the original one.

## API

### Data generation

```@docs
load_test_data
```

### High-level interface
```@docs
IndicatorEvolutionResults
init_metaanalysis_params
analyze_indicators
```

### [Indicators] (@id indicator_functions)
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
threshold_indicators
measure_significance
measure_significances
normalized_confidence_intervall
normalized_percentile_distance
```

### Load data

```@docs
load_test_data()
```