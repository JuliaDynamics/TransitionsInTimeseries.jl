# TransitionIndicators.jl

![TransitionIndicators.jl](https://github.com/JuliaDynamics/JuliaDynamics/blob/master/videos/transitionindicators/logo.gif?raw=true)

```@docs
TransitionIndicators
```

!!! info "Star us on GitHub!"
    If you have found this package useful, please consider starring it on [GitHub](https://github.com/JuliaDynamics/TransitionIndicators.jl).
    This gives us an accurate lower bound of the (satisfied) user count.

## [Content] (@id content)

Multi-stable systems can display abrupt transitions between two stability regimes. To predict such tranistions in real-world systems solely based on data, mathematical tools have been developped in the last decades. Numerous terminologies have been used for them, such as *early warning signals*, *resilience indicators*, *regime-shift identifiers*, *change-point detection* and *transition indicators*. `TransitionIndicators.jl` sticks to the latter terminology and provides an interface that:

- Allows a fast computation of common transition indicators with a couple of lines, as demonstrated in the [example section](@ref example_fastforward).
- Makes the surrogate analysis to test for significance [under the hub](@ref workflow).
- Can be easily extended by any user without touching the source code.
- Reduces the programming overhead for any researcher willing to benchmark new methods.
- Eases the reproducibility thanks to a clear syntax, a simple installation and RNG-seeded surrogate generation.
- Increases trustworthiness thanks to a large test suite.

!!! info "Similar projects"
    An R [toolbox](https://www.early-warning-signals.org/?page_id=42) and a Python [library](https://pypi.org/project/ewstools/) already exist. However, we believe that they are difficult to extend for the user. Furthermore, they do not offer a native performant code, as here allowed by the use of Julia.

## [Approaches] (@id approaches)

Over the last decades, research on transition indicators has largely focused on Critical Slowing Down (CSD). CSD is observed when a system with continuous right-hand side approaches a bifurcation and consists in a resilience loss of the system. For instance this can be diagnosed by an increase of the variance and the AR1-regression coefficient, as demonstrated in the [example section](@ref example_stepbystep). However, we emphasize that this is one out of many possible approaches for obtaining transition indicators. Recent work has explored new approaches relying on nonlinear dynamics or machine learning. `TransitionIndicators.jl` is designed to allow these cutting-edge methods and foster the development of new ones.

## [Under the hub] (@id workflow)

Computing transition indicators is schmetically represented in the plot below and essentially consists of:
1. Detrending the time series to get the fluctuations around the tracked attractor.
1. Estimating the time series of an indicator by sliding a window over the fluctuation time-series.
2. Computing the evolution of the indicator by sliding a window over its time series.
3. Generating many surrogates that preserve important statistical properties of the originial fluctuation time-series.
4. Performing step 2 and 3 for the surrogate time series to check whether the indicator evolution of the original time series shows a significant feature (trend or jump).

The below-depicted fowchart is flexible; boxes can be modified or skipped alltogether!

![Schematic representation of what is happening under the hub.](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/transitionindicators/workflow.svg)

## [Example -- Step-by-step] (@id example_stepbystep)

Let us load data from a bistable nonlinear model subject to noise and to a gradual change of the forcing that leads to a transition. Furthermore, we also load data from a linear model, which is by definition monostable and therefore incapable of transitionning. This is done to control the rate of false positives, a common problem that can emerge when looking for tranisiton indicators. The models are governed by:

```math
\dfrac{\mathrm{d}x_{l}}{\mathrm{d}t} = - x_{l} - 1 + f(t) + n(t) \\
\dfrac{\mathrm{d}x_{nl}}{\mathrm{d}t} = - x_{nl}^3 + x_{nl} + f(t) + n(t)
```

with $x_{l}$ the state of the linear model, $x_{nl}$ the state of the bistable model, $f$ the forcing and $n$ the noise. For $f=0$ they both display an equilibrium point at $x=-1$. However, the bistable model also displays a further equilibrium point at $x=1$. Loading (and visualizing with [`Makie.jl`](https://docs.makie.org/stable/)) such prototypical data to test some indicators can be done by simply running:

```@example MAIN
using TransitionIndicators
using CairoMakie

t, x_linear, x_nlinear = load_linear_vs_doublewell()
X = [x_linear, x_nlinear]
fig, ax = lines(t, x_linear)
lines!(ax, t, x_nlinear)
fig
```

The nonlinear system clearly displays a transition between two stability regimes. To forecast such transition, we analyze the fluctuations of the time series around the tracked attractor. Therefore, a detrending step is needed - here simply obtained by building the difference of the time series with lag 1.

```@example MAIN
fluctuations = [diff(x) for x in [x_linear, x_nlinear]]
t_fluctuations = t[2:end]

fig, ax = lines(t_fluctuations, fluctuations[1])
lines!(ax, t_fluctuations, fluctuations[2])
fig
```

!!! info "Detrending in Julia"
    Detrending can be performed in many ways and therefore remains an external step of `TransitionIndicators.jl`. A wide range of `Julia` packages exists to perform smoothing such as [`Loess.jl`](https://github.com/JuliaStats/Loess.jl) or [`DSP.jl`](https://docs.juliadsp.org/latest/contents/). The detrending step then simply consists of subtracting the smoothed signal from the original one.

We can then compute the time series of indicator by applying a sliding window, determined by the width and the stride with which it is applied. Here we demonstrate this computation with the AR1-regression coefficient (under white-noise assumption), implemented as [`ar1_whitenoise`](@ref):

```@example MAIN
slidingwindow_width = 400
slidingwindow_stride = 1

t_indicator = windowmap(mean, t_fluctuations,
    slidingwindow_width, slidingwindow_stride)
x_indicator_l = windowmap(ar1_whitenoise, fluctuations[1],
    slidingwindow_width, slidingwindow_stride)
x_indicator_nl = windowmap(ar1_whitenoise, fluctuations[2],
    slidingwindow_width, slidingwindow_stride)

fig, ax = lines(t_indicator, x_indicator_l)
lines!(ax, t_indicator, x_indicator_nl)
fig
```

We can obtain the evolution of the indicator over time by applying, here again, a sliding window. We demonstrate this with the computation of the ridge regression slope, a way of measuring a potential increase of the AR1-regression coefficient:

```@example MAIN
slidingwindow_width = 20
slidingwindow_stride = 1
rr = RidgeRegression(t_indicator, slidingwindow_width)

t_evolution = windowmap(mean, t_indicator, slidingwindow_width, slidingwindow_stride)
x_evolution_l = windowmap(rr, x_indicator_l, slidingwindow_width, slidingwindow_stride)
x_evolution_nl = windowmap(rr, x_indicator_nl, slidingwindow_width, slidingwindow_stride)

fig, ax = lines(t_evolution, x_evolution_l)
lines!(ax, t_evolution, x_evolution_nl)
fig
```

As expected from [CSD](@ref approaches), an increase of the AR1-regression coefficient can be observed. Although eyeballing the time series might already be suggestive, surrogates of the fluctuation time series allow to rigorously test the observed increase in AR1-regression coefficient for statistical significance. Furthermore, it allows an automation of the significance computation. Surrogates of the flutuations (here demonstrated for the data generated by the nonlinear model) can be simply generated by:

```@example MAIN
sgen = surrogenerator(fluctuations[2], RandomFourier())
s = sgen()
fig_s, ax_s = lines(t_fluctuations, s)
lines!(ax_s, t_fluctuations, fluctuations[2])
fig_s
```

We can now perform the same analysis for the surrogates. To simplify the syntax, we use [`SignificanceHyperParams`](@ref), a convenience constructor that stores hyperparameters of the significance analysis, such as the sliding-window parameters, the number of surrogaes... etc. Furthermore, it provides default choices for unspecified hyperparameters. We then estimate the trend of the indicators computed on the surrogates by looping over them and using the convenience function [`indicator_evolution`](@ref) that wraps the steps of computing the indicator time series and its evolution. To visualize significant trends, we plot the bands giving the $(-\sigma, \sigma)$, $(-2 \, \sigma, 2 \, \sigma)$ and $(-3 \, \sigma, 3 \, \sigma)$ intervals, with $\sigma$ the standard-deviation of indicator slope across the surrogate time series:

```@example MAIN
p = SignificanceHyperParams(
    n_surrogates = 10_000,      # number of surrogates to test significance
    wv_indicator_width = 400,   # sliding window width for indicator estimation
    wv_indicator_stride = 1,    # sliding window stride for indicator estimation
    wv_evolution_width = 20,    # sliding window width for evolution metric computation
    wv_evolution_stride = 1,    # sliding window stride for evolution metric computation
)

S_evolution = fill(0.0, p.n_surrogates, length(t_evolution))
for i in 1:p.n_surrogates
    s = sgen()
    s_indicator, s_evolution = indicator_evolution(s, ar1_whitenoise, rr, p)
    S_evolution[i, :] .= s_evolution
end
mu = vec(mean(S_evolution, dims = 1))
sigma = vec(std(S_evolution, dims = 1))
band!(ax, t_evolution, mu .- sigma, mu .+ sigma, color = (:red, 0.4) )
band!(ax, t_evolution, mu .- 2 .* sigma, mu .+ 2 .* sigma, color = (:red, 0.2) )
band!(ax, t_evolution, mu .- 3 .* sigma, mu .+ 3 .* sigma, color = (:red, 0.1) )

fig
```

As expected, the data generated by the nonlinear model displays a significant increase of the AR1-regression coefficient before the transition, while the data generated by the linear model does not.

As here demonstrated, performing the step-by-step analysis of transition indicators is possible and might be wishfull for users wanting high flexibility. However, this results in a substantial amount of code. We therefore provide convenience functions that wrap this analysis, as shown in the next section.

## [Example -- Fast-forward] (@id example_fastforward)

!!! tip "Adaptability of `TransitionIndicators.jl`"
    In most cases, computing **various** indicators and evolution metrics are desired. By specifying those either as `Vector{Function}`, the results of this analysis can be computed in a single line with [`analyze_indicators`](@ref). Note that **any** function complying with `f(x::Vector) ➡ y::Real` as an input-output structure can be used here! Most common indicators are already implemented in `TransitionIndicators.jl` and listed [here](@ref indicator_functions) but you can also use user-defined functions! The same holds for the evolution and significance metrics.

The significance analysis against the surrogate time series can subsequently be performed by applying [`measure_significance`](@ref). Finally, if a binary output is desired, one can apply [`threshold_indicators`](@ref):

```@example MAIN
# Initialize the indicator analysis
indicators = [variance, ar1_whitenoise]
evolution_metrics = rr

# Initialize figure and axes
X = [x_linear, x_nlinear]
fig = Figure(resolution = (900, 1100))
nrows, ncols = 5, 2
axs = [Axis(fig[i,j],
    xticklabelsvisible=(i==5 ? true : false)) for i in 1:nrows, j in 1:ncols]
[xlims!(axs[i, j], extrema(t)) for i in 1:nrows, j in 1:ncols]

# Choose which indicator (evolution) to plot in row 3 (and 4)
# Here we choose the first one provided in `indicators`
ind_idx = 1

for j in eachindex(X)
    # Indicator analysis
    result = analyze_indicators(t_fluctuations, fluctuations[j],
        indicators, evolution_metrics, p)
    significance = measure_significance(result, confidence_interval)
    t_indicator = threshold_indicators(result.t_evolution, significance)

    # Plot results
    lines!(axs[1, j], t, X[j])
    lines!(axs[2, j], t[2:end], fluctuations[j])
    lines!(axs[3, j], result.t_indicator, result.X_indicator[1, :,ind_idx])
    lines!(axs[4, j], result.t_evolution, result.X_evolution[1, :, ind_idx])
    lines!(axs[5, j], result.t_evolution, significance[1, :, 1])
    lines!(axs[5, j], result.t_evolution, significance[1, :, 2])
    [vlines!(axs[i, j], t_indicator, linestyle = :dash, color = :red) for i in [1, 5]]
end
fig
```
Here we see that a transition is correctly forecasted for the bistable system, whereas none is predicted for the linear system.

!!! warning "Thresholding significance"
    Although thresholding the output of a significance computation might be hard to avoid in some applications, we recommend to rather look at the significance time series, as they provide richer, non-binary information.

## API

### Load data

```@docs
load_linear_vs_doublewell()
```

### High-level interface
```@docs
IndicatorEvolutionResults
SignificanceHyperParams
analyze_indicators
indicator_evolution
windowmap
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
RidgeRegression
kendalltau
spearman
```

### Surrogates

For the surrogate generation, you can use any subtype of `Surrogate` defined in [TimeseriesSurrogates.jl](https://juliadynamics.github.io/TimeseriesSurrogates.jl/stable/#Surrogate-methods).

### Significance metrics

```@docs
threshold_indicators
measure_significance
confidence_interval
normalized_percentile
```