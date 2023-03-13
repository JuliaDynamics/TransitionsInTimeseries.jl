
## [Tutorial -- Educational] (@id example_stepbystep)

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
indicator = ar1_whitenoise
indicator_window = (width = 400, stride = 1)

t_indicator = windowmap(midpoint, tfluct; indicator_window...)

indicator_l = windowmap(indicator, x_l_fluct; indicator_window...)

indicator_nl = windowmap(indicator, x_nl_fluct; indicator_window...)

fig, ax = lines(t_indicator, indicator_l)
lines!(ax, t_indicator, indicator_nl)
ax.title = "indicator timeseries"
fig
```

The lines plotted above are the **indicator timeseries**.

### Change metric timeseries

From here, we process the **indicator timeseries** to quantify changes in it. This step is in essence the same as before: we apply some function over a sliding window of the indicator timeseries. We call this new timeseries the **change metric timeseries**. In the example here, the change metric we will employ will be the slope (over a sliding window), calculated via means of a ridge regression


```@example MAIN
change_window = (width = 20, stride = 1)
ridgereg = RidgeRegression(t_indicator, change_window.width)

t_change = windowmap(midpoint, t_indicator; change_window...)
change_l = windowmap(ridgereg, indicator_l; change_window...)
change_nl = windowmap(ridgereg, indicator_nl; change_window...)

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
indicator_s = windowmap(indicator, s; indicator_window...)
change_s = windowmap(ridgereg, indicator_s; change_window...)

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
        indicator_s = windowmap(indicator, s; indicator_window...)
        change_s = windowmap(ridgereg, indicator_s; change_window...)
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

## [Tutorial -- TransitionIndicators.jl] (@id example_fastforward)

TransitionIndicators.jl wraps the process of the step-by-step example into a simple, extendable, and modular API that researchers can use. In addition, it allows performing the same analysis for several different indicators / change metrics in one go.

The workflow is simple and it works is as follows:

1. Create an instance of [`IndicatorsConfig`](@ref), that dictates which indicators will be used, and over what sliding window.
2. Create an instance of [`SignificanceConfig`](@ref), that dictates what change metrics, surrogate types, and sliding window will be used to quantify a significant change of the indicators.
3. Along with the input timeseries `x` these three are plugged into [`indicators_analysis`](@ref).
4. The output, which is an [`IndicatorsResults`](@ref), can be used to visualize the results, or given to [`indicators_significance`](@ref) to provide flags of the timepoints where there is a significant change of each indicator.

The following blocks apply this process, and visualize it, for the examples we used in the educational part of the tutorial.


First we load, and do the necessary pre-processing, to get input data
```@example MAIN
using TransitionIndicators, CairoMakie

t, x_linear, x_nlinear = load_linear_vs_doublewell()

x_nl_fluct = diff(x_nlinear)
tfluct = t[2:end]

fig, ax = lines(tfluct, x_l_fluct)
ax.title = "input timeseries"
fig
```

Then we decide what indicators and change metrices to use

```@example MAIN
# these indicators are suitable for Critical Slowing Down
indicators = [var, ar1_whitenoise]
ind_conf = IndicatorsConfig(indicators; width = 400)

# use spearman correlation for both indicators
change_metric = spearman
sig_conf = SignificanceConfig(change_metric; width = 20, n_surrogates = 1000)

# perform the full analysis
result = indicators_analysis(x_nl_fluct, ind_conf, sig_conf)
```

And lastly, obtain some flags for when there is a significant indicator change
```@example MAIN
sig = indicators_significance(result, 0.99)

# Plot the original timeseries
fig, ax = lines(tfluct, x_nl_fluct; color = Cycled(2), label = "input")
# Scatter the significance of each indicator
for i in size(sig, 2)
    signif_idxs = findall(sig[:, i])
    isempty(signif_idxs) && continue
    # get timepoints in real time
    tflags = tfluct[result.t_change[signif_idxs]]
    vlines!(ax, tflags; label = "indicator $(indicators[i])")
end
axislegend()
fig
```
