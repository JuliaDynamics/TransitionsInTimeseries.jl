# TransitionIndicators.jl

![TransitionIndicators.jl](https://github.com/JuliaDynamics/JuliaDynamics/blob/master/videos/transitionindicators/logo.gif?raw=true)

```@docs
TransitionIndicators
```

!!! info "Star us on GitHub!"
    If you have found this package useful, please consider staridgereging it on [GitHub](https://github.com/JuliaDynamics/TransitionIndicators.jl).
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

## [Under the hood] (@id workflow)

Computing transition indicators is schmetically represented in the plot below and essentially consists of:
1. Detrending the time series to get the fluctuations around the tracked attractor.
1. Estimating the time series of an indicator by sliding a window over the fluctuation time-series.
2. Computing the evolution of the indicator by sliding a window over its time series.
3. Generating many surrogates that preserve important statistical properties of the originial fluctuation time-series.
4. Performing step 2 and 3 for the surrogate time series to check whether the indicator evolution of the original time series shows a significant feature (trend or jump).

The below-depicted fowchart is flexible; boxes can be modified or skipped alltogether!

![Schematic representation of what is happening under the hub.](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/transitionindicators/workflow.svg)

## [Example -- Step-by-step] (@id example_stepbystep)

### Raw input data

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
X = (x_linear, x_nlinear)
fig, ax = lines(t, x_linear)
lines!(ax, t, x_nlinear)
ax.title = "raw data"
fig
```

### Preprocessing

!!! note "Not part of TransitionIndicators.jl"
    Any timeseries pre-processing, such as the de-trending step we do here,
    is not part of TransitionIndicators.jl and is the responsibility of the researcher.


The nonlinear system clearly displays a transition between two stability regimes. To forecast such transition, we analyze the fluctuations of the time series around the tracked attractor. Therefore, a detrending step is needed - here simply obtained by building the difference of the time series with lag 1.

```@example MAIN
x_l_fluct = diff(x_linear)
x_nl_fluct = diff(x_nlinear)
tfluct = t[2:end]

fig, ax = lines(tfluct, x_l_fluct)
lines!(ax, tfluct, x_nl_fluct .+ 0.05)
ax.title = "input timeseries"
fig
```

At this point, `x_l_fluct` and `x_nl_fluct` are considered the **input timeseries**.

!!! info "Detrending in Julia"
    Detrending can be performed in many ways and therefore remains an external step of `TransitionIndicators.jl`. A wide range of `Julia` packages exists to perform smoothing such as [`Loess.jl`](https://github.com/JuliaStats/Loess.jl) or [`DSP.jl`](https://docs.juliadsp.org/latest/contents/). The detrending step then simply consists of subtracting the smoothed signal from the original one.

### Indicator timeseries

We can then compute the values of some "indicator" (a Julia function that inputs a timeseries and outputs a number). An indicator should be a quantity that is likely to change if a transition occurs in the timeseries. We compute indicators by applying a sliding window over the **input timeseries**, determined by the width and the stride with which it is applied. Here we demonstrate this computation with the AR1-regression coefficient (under white-noise assumption), implemented as [`ar1_whitenoise`](@ref):

```@example MAIN
width_indicator = 400
stride_indicator = 1
indicator = ar1_whitenoise

t_indicator = windowmap(mean, tfluct, width_indicator, stride_indicator)

indicator_l = windowmap(indicator, x_l_fluct, width_indicator, stride_indicator)

indicator_nl = windowmap(indicator, x_nl_fluct, width_indicator, stride_indicator)

fig, ax = lines(t_indicator, indicator_l)
lines!(ax, t_indicator, indicator_nl)
ax.title = "indicator timeseries"
fig
```

The lines plotted above are the **indicator timeseries**.

### Change metric timeseries

From here, we process the **indicator timeseries** to quantify changes in it. This step is in essence the same as before: we apply some function over a sliding window of the indicator timeseries. We call this new timeseries the **change metric timeseries**. In the example here, the change metric we will employ will be the slope (over a sliding window), calculated via means of a ridge regression


```@example MAIN
width_change = 20
stride_change = 1
ridgereg = RidgeRegression(t_indicator, width_change)

t_change = windowmap(mean, t_indicator, width_change, stride_change)
change_l = windowmap(ridgereg, indicator_l, width_change, stride_change)
change_nl = windowmap(ridgereg, indicator_nl, width_change, stride_change)

fig, ax = lines(t_change, change_l)
lines!(ax, t_change, change_nl)
ax.title = "change metric timeseries"
fig
```


### Timeseries surrogates

As expected from [Critical Slowing Down](@ref approaches), an increase of the AR1-regression coefficient can be observed. Although eyeballing the time series might already be suggestive, we want a rigorous framework for testing for significance.

In TransitionsIdentifiers.jl we perform significance testing using the method of timeseries surrogates and the [TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl) Julia package. This has the added benefits of flexibility in choosing the surrogate generation method, reproducibility, and automation. Note that `TimeseriesSurrogates` is re-exported by `TransitionIndicators`, so that you don't have to `using` both of them.

To illustrate the surrogate, we compare the change metric computed from the bistable timeseries what that computed from a surrogate of the same timeseries.

```@example MAIN
# Generate Fourier random-phase surrogates
using Random: Xoshiro
s = surrogate(x_nl_fluct, RandomFourier(), Xoshiro(123))
fig, ax = lines(tfluct, x_nl_fluct; color = Cycled(2))
lines!(ax, tfluct, s .- 0.05; color = Cycled(3))
ax.title = "real signal vs. surrogate(s)"

# compute and plot change metric
indicator_s = windowmap(indicator, s, width_indicator, stride_indicator)
change_s = windowmap(ridgereg, indicator_s, width_change, stride_change)

ax, = lines(fig[1,2], t_change, change_nl; color = Cycled(2), label = "nonlin")
lines!(ax, t_change, change_s; color = Cycled(3), label = "surrogate")
axislegend()
ax.title = "change metric"

fig
```

### Quantifying significance

To quantify the significance of the values of the **change metric timeseries** we perform a standard surrogate test. We calculate the change metric for thousands of surrogates of the input timeseries, and then detect the points in time where the change metric timeseries exceeds a threshold (such as the 95-quantile) of the change metrics of the surrogates.

To visualize significant trends, we plot the bands giving the $(-2 \, \sigma, 2 \, \sigma)$ and $(-3 \, \sigma, 3 \, \sigma)$ intervals of the surrogate values, with $\sigma$ the standard-deviation of indicator slope across the surrogate time series:

```@example MAIN
n_surrogates = 1_000
fig = Figure()
axl = Axis(fig[1,1]; title = "linear")
axnl = Axis(fig[1,2]; title = "nonlinear")

for (j, ax, x) in zip(1:2, (axl, axnl), (x_l_fluct, x_nl_fluct))

    sgen = surrogenerator(x, RandomFourier(), Xoshiro(123))
    change_s_distr =  zeros(n_surrogates, length(change_s))

    # Collect all surrogate change metrics
    for i in 1:n_surrogates
        s = sgen()
        indicator_s = windowmap(indicator, s, width_indicator, stride_indicator)
        change_s = windowmap(ridgereg, indicator_s, width_change, stride_change)
        change_s_distr[i, :] .= change_s
    end

    mu = vec(mean(change_s_distr, dims = 1))
    sigma = vec(std(change_s_distr, dims = 1))

    # Plot (real signal) change metric and various confidence intervals
    orig_change = j == 1 ? change_l : change_nl
    lines!(ax, t_change, orig_change; color = Cycled(j))
    band!(ax, t_change, mu .- 2 .* sigma, mu .+ 2 .* sigma, color = (:red, 0.2) )
    band!(ax, t_change, mu .- 3 .* sigma, mu .+ 3 .* sigma, color = (:red, 0.1) )
end

fig
```

As expected, the data generated by the nonlinear model displays a significant increase of the AR1-regression coefficient before the transition, while the data generated by the linear model does not.

Performing the step-by-step analysis of transition indicators is possible and might be preferred for users wanting high flexibility. However, this results in a substantial amount of code. We therefore provide convenience functions that wrap this analysis, as shown in the next section.

## [Example -- Fast-forward] (@id example_fastforward)

!!! tip "Adaptability of `TransitionIndicators.jl`"
    In most cases, computing **various** indicators and evolution metrics are desired. By specifying those either as `Vector{Function}`, the results of this analysis can be computed in a single line with [`analyze_indicators`](@ref). Note that **any** function complying with `f(x::Vector) → y::Real` as an input-output structure can be used here! Most common indicators are already implemented in `TransitionIndicators.jl` and listed [here](@ref indicator_functions) but you can also use user-defined functions! The same holds for the evolution and significance metrics.

The significance analysis against the surrogate time series can subsequently be performed by applying [`measure_significance`](@ref). Finally, if a binary output is desired, one can apply [`threshold_indicators`](@ref):

```@example MAIN
# Initialize the indicator analysis
indicators = [variance, ar1_whitenoise]
evolution_metrics = ridgereg

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
    t_indicator = threshold_indicators(result.t_change, significance)

    # Plot results
    lines!(axs[1, j], t, X[j])
    lines!(axs[2, j], t[2:end], fluctuations[j])
    lines!(axs[3, j], result.t_indicator, result.X_indicator[1, :,ind_idx])
    lines!(axs[4, j], result.t_change, result.X_evolution[1, :, ind_idx])
    lines!(axs[5, j], result.t_change, significance[1, :, 1])
    lines!(axs[5, j], result.t_change, significance[1, :, 2])
    [vlines!(axs[i, j], t_indicator, linestyle = :dash, color = :red) for i in [1, 5]]
end
fig
```
Here we see that a transition is coridgeregectly forecasted for the bistable system, whereas none is predicted for the linear system.

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

### surrogates

For the surrogate generation, you can use any subtype of `surrogate` defined in [Timeseriessurrogates.jl](https://juliadynamics.github.io/Timeseriessurrogates.jl/stable/#surrogate-methods).

### Significance metrics

```@docs
threshold_indicators
measure_significance
confidence_interval
normalized_percentile
```