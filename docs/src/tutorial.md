# Tutorial

## [Workflow] (@id workflow)

Computing transition indicators consists of the following steps:

1. Doing any pre-processing of raw data first such as detrending (_not part of TransitionIndicators.jl_). This yields the **input timeseries**.
2. Estimating the timeseries of an indicator by sliding a window over the input timeseries.
3. Computing the changes of the indicator by sliding a window over its timeseries.
4. Generating many surrogates that preserve important statistical properties of the original timeseries.
5. Performing step 2 and 3 for the surrogate timeseries.
6. Checking whether the indicator change timeseries of the real timeseries shows a significant feature (trend, jump or anything else) when compared to the surrogate data.

These steps are illustrated one by one in the tutorial below, and then summarized in the convenient API that TransitionIndicators.jl exports.

## [Tutorial -- Educational] (@id example_stepbystep)

### Raw input data

Let us load data from a bistable nonlinear model subject to noise and to a gradual change of the forcing that leads to a transition. Furthermore, we also load data from a linear model, which is by definition monostable and therefore incapable of transitioning. This is done to control the rate of false positives, a common problem that can emerge when looking for transition indicators. The models are governed by:

```math
\dfrac{\mathrm{d}x_{l}}{\mathrm{d}t} = - x_{l} - 1 + f(t) + n(t) \\
\dfrac{\mathrm{d}x_{nl}}{\mathrm{d}t} = - x_{nl}^3 + x_{nl} + f(t) + n(t)
```

with $x_{l}$ the state of the linear model, $x_{nl}$ the state of the bistable model, $f$ the forcing and $n$ the noise. For $f=0$ they both display an equilibrium point at $x=-1$. However, the bistable model also displays a further equilibrium point at $x=1$. Loading (and visualizing with [Makie](https://docs.makie.org/stable/)) such prototypical data to test some indicators can be done by simply running:

```@example MAIN
using TransitionIndicators, CairoMakie

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


The nonlinear system clearly displays a transition between two stability regimes. To forecast such transition, we analyze the fluctuations of the timeseries around the tracked attractor. Therefore, a detrending step is needed - here simply obtained by building the difference of the timeseries with lag 1.

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
    Detrending can be performed in many ways. A wide range of Julia packages exists to perform smoothing such as [Loess.jl](https://github.com/JuliaStats/Loess.jl) or [DSP.jl](https://docs.juliadsp.org/latest/contents/). There the detrending step consists of subtracting the smoothed signal from the original one.

### Indicator timeseries

We can then compute the values of some "indicator" (a Julia function that inputs a timeseries and outputs a number). An indicator should be a quantity that is likely to change if a transition occurs in the timeseries. We compute indicators by applying a sliding window over the **input timeseries**, determined by the width and the stride with which it is applied. Here we demonstrate this computation with the AR1-regression coefficient (under white-noise assumption), implemented as [`ar1_whitenoise`](@ref):

```@example MAIN
indicator = ar1_whitenoise
indicator_window = (width = 400, stride = 1)

# left-bracketing respects information available until t. Alternatives: center, right.
t_indicator = slidebracket(tfluct, :left; indicator_window...)
indicator_l = windowmap(indicator, x_l_fluct; indicator_window...)
indicator_nl = windowmap(indicator, x_nl_fluct; indicator_window...)

fig, ax = lines(t_indicator, indicator_l)
lines!(ax, t_indicator, indicator_nl)
ax.title = "indicator timeseries"
fig
```

The lines plotted above are the **indicator timeseries**.

### Change metric timeseries

From here, we process the **indicator timeseries** to quantify changes in it. This step is in essence the same as before: we apply some function over a sliding window of the indicator timeseries. We call this new timeseries the **change metric timeseries**. In the example here, the change metric we will employ will be the slope (over a sliding window), calculated via means of a [`RidgeRegressionSlope`](@ref):


```@example MAIN
change_window = (width = 20, stride = 1)
cmp = ChangeMetricsParams(lambda_ridge = 0.0)   # can be set >0 to have regularization
ridgereg = RidgeRegressionSlope(t_indicator[1:change_window.width], cmp)

t_change = slidebracket(t_indicator, :left; change_window...)
change_l = windowmap(ridgereg, indicator_l; change_window...)
change_nl = windowmap(ridgereg, indicator_nl; change_window...)

fig, ax = lines(t_change, change_l)
lines!(ax, t_change, change_nl)
ax.title = "change metric timeseries"
fig
```

### Timeseries surrogates

As expected from [Critical Slowing Down](@ref approaches), an increase of the AR1-regression coefficient can be observed. Although eyeballing the timeseries might already be suggestive, we want a rigorous framework for testing for significance.

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

To quantify the significance of the values of the **change metric timeseries** we perform a standard surrogate test by computing the [p-value](https://en.wikipedia.org/wiki/P-value) w.r.t. the change metrics of thousands of surrogates of the input timeseries. A low p-value (typically `p<0.05`) is commonly considered as significant. To visualize significant trends, we plot the p-value vs. time:

```@example MAIN
n_surrogates = 1_000
fig = Figure()
axl = Axis(fig[1,1]; title = "linear")
axnl = Axis(fig[1,2]; title = "nonlinear")
axsigl = Axis(fig[2,1])
axsignl = Axis(fig[2,2])

for (j, ax, axsig, x) in zip(1:2, (axl, axnl), (axsigl, axsignl), (x_l_fluct, x_nl_fluct))

    orig_change = j == 1 ? change_l : change_nl
    sgen = surrogenerator(x, RandomFourier(), Xoshiro(123))
    pval = zeros(length(change_s))

    # Collect all surrogate change metrics
    for i in 1:n_surrogates
        s = sgen()
        indicator_s = windowmap(indicator, s; indicator_window...)
        change_s = windowmap(ridgereg, indicator_s; change_window...)
        pval += orig_change .< change_s
    end

    pval ./= n_surrogates
    lines!(ax, t_change, orig_change; color = Cycled(j))
    lines!(axsig, t_change, pval; color = Cycled(j+2))
end

fig
```

As expected, the data generated by the nonlinear model displays a significant increase of the AR1-regression coefficient before the transition, which is manifested by a low p-value. In contrast, the data generated by the linear model does not show anything similar.

Performing the step-by-step analysis of transition indicators is possible and might be preferred for users wanting high flexibility. However, this results in a substantial amount of code. We therefore provide convenience functions that wrap this analysis, as shown in the next section.

## [Tutorial -- TransitionIndicators.jl] (@id example_fastforward)

TransitionIndicators.jl wraps this typical workflow into a simple, extendable, and modular API that researchers can use with little effort. In addition, it allows performing the same analysis for several indicators / change metrics in one go.

The interface is simple, and directly parallelizes the [Workflow](@ref) as follows:

1. Create an instance of [`IndicatorsConfig`](@ref), that dictates which indicators will be used, and over what sliding window.
2. Create an instance of [`SignificanceConfig`](@ref), that dictates what change metrics, surrogate types, and sliding window will be used to quantify a significant change of the indicators.
3. Along with the input timeseries `x` these three are plugged into [`indicators_analysis`](@ref).
4. The output, which is an [`IndicatorsResults`](@ref), can be used to visualize the results.

The following blocks apply this process, and visualize it, for the examples we used in the educational part of the tutorial. First we load, and do the necessary pre-processing, to get input data:

```@example MAIN
using TransitionIndicators, CairoMakie

t, x_linear, x_nlinear = load_linear_vs_doublewell()

x_nl_fluct = diff(x_nlinear)
tfluct = t[2:end]

fig, ax = lines(tfluct, x_nl_fluct)
ax.title = "input timeseries"
fig
```

Then we decide what indicators and change metrics to use and run the analysis.

```@example MAIN
# these indicators are suitable for Critical Slowing Down
indicators = [var, ar1_whitenoise]
indconfig = IndicatorsConfig(tfluct, indicators; width = 400)

# use the ridge regression slope for both indicators
change_metrics = [RidgeRegressionSlope]
sigconfig = SignificanceConfig(indconfig, change_metrics; width = 30, n_surrogates = 1000)

# perform the full analysis
results = indicators_analysis(t, x_nl_fluct, indconfig, sigconfig)
```

That's it, just a couple of lines! To get insight on the results, we plot the obtained p-values vs. the original time series:

```@example MAIN
fig, ax = lines(tfluct, x_nl_fluct; color = Cycled(2), label = "input")
axpval, = lines(fig[2,1], sigconfig.t_change, results.pval[:, 1]; color = Cycled(3), label = "p-value of var")
lines!(axpval, sigconfig.t_change, results.pval[:, 2]; color = Cycled(4), label = "p-value of ar1")
xlims!(ax, (0, 50))
xlims!(axpval, (0, 50))
axislegend(axpval)
fig
```

Thresholding the p-value results in a loss of information and is therefore to be avoided. For automation purposes it might however be unavoidable. Here we show an example where we consider that both the variance and the AR1-regression coefficient need to show a p-value below 0.05:

```@example MAIN
threshold = 0.05
signif_idxs = vec(count(results.pval .< threshold, dims = 2) .>= 2)
tflags = sigconfig.t_change[signif_idxs]
vlines!(ax, tflags; label = "flags", color = Cycled(3), linestyle = :dash)
axislegend(ax)
fig
```