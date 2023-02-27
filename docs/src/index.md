# TransitionIndicators.jl

```@docs
TransitionIndicators
```

## Why do we need TransitionIndicators.jl?

Multi-stability is an ubiquitous feature of real-world dynamical systems, as for instance in epidemiology, sociology and climate science. The transition between two stability regimes can be abrupt and of large magnitude, although its trigger can be gradual and of small amplitude. The repercusions of such transitions can be hard to mitigate and finding mathematical tools to predict their occurence solely based on data has therefore been an active field of research over the last decades. Hereby, many terms were employed such as "early warning signals", "resilience indicators", "regime-shift identifiers", "change-point detection" and "transition indicators". The present package will use the latter terminology.

As research on transition indicators gains interest, many problems are encountered:
- applying state-of-the-art techniques may be computationally expensive, especially when studying data in time and space
- developping new indicators requires researchers to programm some standard functionalities, thus unnecessarily increasing the overhead of their investigation
- benchmarking a new transition indicator requires the implementation of pre-existing techniques, thus further increasing the aformentionned overhead
- some published results are hard to reproduce

Worth mentioning: an `R` toolbox and a `Python` library addressing some of these problems already exist. However, we believe that they are difficult to extend for the user and they do not offer a native performant code, as allowed by the use of `julia`.

TransitionIndicators.jl aims to tackle all of these problems by providing an optimized code with a high-level interface that can be easily extended for studying new transition indicators. It comes  with most of the existing transition indicators already implemented to ease benchmarking. We hope it to contribute to reproducibility by its simple installation:

```julia
] add TransitionIndicators
using TransitionIndicators
```

More than anything, `TransitionIndicators.jl` allows to perform a full analysis of transition indicators in few lines, which is demonstrated in the  section.

## What is TransitionIndicators.jl computing exaclty?



## [Example] 

Let us generate data from a bi-stable nonlinear model subject to noise and to a gradual change of the forcing that leads to a transition. Furthermore, we also generate data from a linear model, which is by definition monostable and therefore incapable of transitionning. This is done to control the rate of false positives. Generating (and plotting) such prototypical data to test some indicators can be done by simply running:

```@example MAIN
using TransitionIndicators
using CairoMakie

t, x_linear, x_nlinear = generate_test_data()
fig = Figure()
axs = [Axis(fig[i, j]) for i in 1:5, j in 1:2]
lines!(axs[1, 1], t, x_linear)
lines!(axs[1, 2], t, x_nlinear)
fig
```

Finding transition indicators focuses on analyzing the fluctuations of the system around the tracked attractor. Therefore a detrending step is common and can be performed by the user, for instance with the help of a wide range of `julia` packages such as `DSP.jl`, `SmoothingSplines.jl` or `Loess.jl`. Here we use the latter one:

```
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

```
p = init_metaanlysis_params()

indicators = [var, ar1_whitenoise]
m = precompute_ridge_slope(1.0:p.wv_evolution_width+1, lambda = 0.1)
evolution_metrics = curry(precomputed_ridge_slope, m)

results_lin = analyze_indicators(t, lin_fluctuations, indicators, evolution_metrics, p)
results_nonlin = analyze_indicators(t, nlin_fluctuations, indicators, evolution_metrics, p)
```

### Copy-pastable code

```
model = loess(xs, ys, span=0.5)
```


using Statistics: mean, var

n = 1001
t = collect(1.0:n)
x = copy(t)
p = init_metaanalysis_params(n_surrogates = 100)
res = analyze_indicators(t, x, [mean, var], ridge_slope, p)

# The trend of mean(windowview) is the stride for x=t
meantrend_ground_truth = fill(p.wv_indicator_stride, length(res.t_evolution))
# The trend of var(windowview) is 0 for x any affine function of t.
vartrend_ground_truth = fill(0.0, length(res.t_evolution))

## API

```@docs
analyze_indicators
```
