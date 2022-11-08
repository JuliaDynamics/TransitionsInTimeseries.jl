# [Example: double-well vs. linear system](@id doublewell_vs_linear)

## Introduction

One of the simplest system displaying bistability is the double-well potential. By taking the negative gradient of such 1D potential, one obtains following normal form:

```math
\dot{x} = -x^3 + x + F(t) + n(t)
```

with $\dot{x}$ the time derivative of the state $x$, $F$ the forcing and $n$ a source of additive white noise with variance $\sigma^2$. For the autonomous case, the system displays stable equilibrium points in $x=-1$, $x=1$ and an unstable one in $x=0$. By setting $x(t=0) = -1$ and applying a slow linear drift in $F$, one can artificially generate a transition, which is expected to be detectable thanks to TIs.

A common problem of TIs is to deliver false positives, i.e. prediction of a transition ahead even when none is about to happen. To investigate this potential drawback in parallel, we study a linear model also displaying an equilibrium point at $x=-1$ for the autonomous case and defined as:

```math
\dot{x} = \lambda (x+1) + F(t) + n(t).
```

Here we fix $\lambda=-1$. As linear models do not allow critical transitions to occur, they provide a baseline to check for false positives. Putting these equations into code results in:

```julia
T = Float32

# Deterministic part of linear system.
function f_linear(dx, x, p, t)
    dx[1] = p["λ"] * (x[1] + 1) + forcing(p, t)
end

# Deterministic part of double-well system.
function f_doublewell(dx, x, p, t)
    dx[1] = (- x[1]^3 + x[1]) + forcing(p, t)
end

# Stochastic part of both system for noise.
function g_whitenoise(dx, x, p, t)
    dx[1] = p["σ"]
end

# Linear forcing
function forcing(p, t)
    return p["α"] * t
end

pmodel = Dict(
    "σ" => T(0.02),
    "λ" => T(-1),
    "α" => T(0.01),
)
```

## Generating data of critical transition

Let's initialize both system at their common equilibrium point $x=-1$ and generate forced time series by integrating the resulting stochastic differential equation.

```julia
using DifferentialEquations

dt = 5f-2
t = collect(0f0:dt:100f0)
tspan = extrema(t)
x0 = [-1f0]

models = [f_linear, f_doublewell]
labels = string.(models)
nx = length(models)
X = zeros(T, nx, length(t))

for (f, i) in zip(models, 1:nx)
    prob = SDEProblem(f, g_whitenoise, x0, tspan, pmodel)
    sol = solve(prob, EM(), dt=dt)
    X[i, :] = vcat(sol.u...)
end
```

In the last step, we integrated both problems by an Euler-Maruyama scheme provided in [`DifferentialEquations.jl`](https://diffeq.sciml.ai/stable/tutorials/sde_example/). The solutions were stored in a matrix $X$ with following structure:

```math
X = \begin{pmatrix}
& x^{(1)}_1 & x^{(1)}_2 & \dots & x^{(1)}_{n_t- 1} & x^{(1)}_{n_t} & \\
& \\
& x^{(2)}_1 & x^{(2)}_2 & \dots & x^{(2)}_{n_t- 1} & x^{(2)}_{n_t} &
\end{pmatrix},
```

where $x^{(j)}_k$ the state of model $j$ (1=linear, 2=double-well) at time step $k$. Let us briefly visualise the time series:

```julia
using CairoMakie

ylabels = [
    L"$x(t)$",
    L"residual $r(t)$",
    L"AR1 regression coefficient $\hat{\theta}(t)$",
    L"Trend $\alpha(t)$ of $\hat{\theta}(t)$",
    L"Percentile significance $p$ of $\alpha(t)$",
]

nrows, ncols = length(ylabels), nx
fig = Figure( resolution = (1500, 1500) )
axs = [[Axis(
    fig[i,j],
    title = ( i==1 ? labels[j] : " " ),
    xlabel = ( i==nrows ? "Time" : " " ),
    ylabel = ( j==1 ? ylabels[i] : " " ),
    xminorticks = tspan[1]:5:tspan[2],
    xminorgridvisible = true,
    ) for j in 1:ncols] for i in 1:nrows]

for j in 1:nx
    lines!(axs[1][j], t, X[j, :], label = L"data $\,$")
    ylims!(axs[1][j], (-1.2, 1.5))
end
fig
```

## Smoothing and de-trending

Critical slowing down is possibly observed in the response of the system to noise. To study this without taking the slower, deterministic dynamic into account, the time series is filtered. Here some of the many options are:

1. Take a moving average.

2. Perform a kernel smoothing (convolution of the signal with a kernel). Note that the previously mentioned moving average is a smoothing with uniform kernel and is therefore a special case of kernel smoothing.

3. Apply a difference operation.

4. Apply a low-pass filter to get the trend and remove it from the original time series. 

5. Apply a high-pass filter to directly get the faster dynamic. Whereas this method is more direct than the previous one, it might be harder to monitor the level of smoothing.

Whereas kernel smoothing and difference operations are integrated to the present package, filter design and application is to be performed by the user, for instance by using [`DSP.jl`](https://docs.juliadsp.org/stable/filters/).

In the case of a smoothing with discrete kernel $u \in \mathbb{R}^{2 N+1}$, a sliding window is applied on the time series $x$ to compute the smoothed values $\tilde{x}$, which can be expressed as following dot product:

```math
\tilde{x}_{i} = \begin{pmatrix} & x_{i-N} & x_{i-N+1} & \dots &
 x_{i+N-1} & x_{i+N} & \end{pmatrix} \cdot u
```

An extensive list of kernels can be found [here](https://en.wikipedia.org/wiki/Kernel_(statistics)). They are implemented in the present package and can be found in the [API reference](@ref api_ref). In the present case, a centered window is applied on the data $x$, but the user can choose whether they prefer to use e.g. a left-window. Furthermore, the half-width $N$ of the window as well as the stride with which the operation should be applied are to be determined by the user. An exemplary code is here provided:

```julia
using TransitionIdentifiers

# half-width and stride of smoothing window.
T_smooth_wndw = T(10 * dt)
T_smooth_strd = dt
p_window_smooth = get_windowing_params([dt, T_smooth_wndw, T_smooth_strd])

window = centered_wndw  # alternative: left_window, right_window
kernel = uniform_kernel # alternative: gaussian_kernel and many more.

Xtrend = get_trend( X, p_window_smooth, window, kernel )
ttrend = trim_wndw( t, p_window_smooth, window )

# Compute the residuals of the trend fitting.
Xres = trim_wndw(X, p_window_smooth, window) - Xtrend

for i in 1:nx
    lines!(axs[1][i], ttrend, Xtrend[i, :], label = L"trend $\,$")
    lines!(axs[2][i], ttrend, Xres[i, :], label = L"original residual $\,$")
end

axislegend(axs[1][2], position = :lt)
fig
```

The result of the detrended time series is called _residual_.

## A step-by-step introduction

For didactic purposes we now show step by step how to compute a transition indicator (here exemplarly the AR1 regression coefficient under white noise assumption), its trend and its significance. This way of proceeding is more verbose in syntax but lets the user have full control of possible in-between steps they might want to perform. A faster syntax wrapping all these steps is shown in [The fast-forward way](@ref ffway).

### Estimating the transition indicator

Now that the residuals are computed, one can estimate transition indicators. For this we first need to define windowing parameters and compute the associated time vector. The computation itself is performed via the [`slide_estimator`](@ref) function, which slides the windowed computation over the time series.

```julia
# Windowing parameters for TI estimation.
T_indctr_wndw = T(5e-1)
T_indctr_strd = T(dt)
p_window_indctr = get_windowing_params([dt, T_indctr_wndw, T_indctr_strd])
window = left_wndw

tindctr = trim_wndw(ttrend, p_window_indctr, window)
origin_ar1 = slide_estimator(Xres, p_window_indctr, ar1_whitenoise, window)

for j in 1:nx
    lines!(axs[3][j], tindctr, origin_ar1[j, :], label = L"AR1 coefficient $\,$")
    ylims!(axs[3][j], (-.1, 1.1))
end

fig
```

The implemented functions to estimate other transition indicators are listed in the [API reference](@ref api_indicators). They can be used instead of `ar1_whitenoise` while keeping the rest of the syntax identical.

### Estimating indicator trends

To indicate a transition, one is primarly interested in the trend of an indicator rather than in its absolute value. Here again, the computation requires the definition of a certain window over which the trend computation is performed. To estimate the trend of a variable, an affine regression can be performed. The coefficient of the linear term gives a scalar measure of the trend. The ridge regression is a generalisation of the affine regression as it allows for regularization:

```julia
T_idtrend_wndw = T(2e0)
T_idtrend_strd = T(1e0)
dt_id = round(tindctr[2] - tindctr[1]; digits = 4)
p_window_idtrend = get_windowing_params([dt_id, T_idtrend_wndw, T_idtrend_strd])

tidtrend = trim_wndw( tindctr, p_window_idtrend, window )
origin_ar1_trend = slide( origin_ar1, tindctr, p_window_idtrend, ridge_regression_slope, window)

for j in 1:nx
    lines!(axs[4][j], tidtrend, origin_ar1_trend[j, :], label = L"trend of AR1 coefficient $\,$")
end
fig
```

Many more methods are available to estimate the trend of a time series. They are listed in the [API reference](@ref api_trends) and can be used instead of `ridge_regression_slope` while using the same syntax.

### Significance of indicator trends

It is of course possible to plot the trend over time and eyeball any large value to recognise a possible transition. However this presents several major drawbacks:

- Such arbitrary criterion is not well-defined.
- It does not allow automation, which is highly desired when dealing with large data.

A solution to this problem is to perform a test for statistical significance. To this end, we artificially generate so-called surrogate time series (typically 1,000 -- 10,000 for each original time series). These preserve most properties of the original residual, while suppressing its deterministic content. For a given indicator, its increase represents a certain percentile of the surrogate data, if this percentile is high enough, the indicator increase can be considered as significant.

Generating surrogates can be done by running a single command. Furthermore, we here visualise some of the surrogate residuals to give the user following intuition: although the spectrum is not altered compared to the original residuals, Fourier surrogates do not display any deterministic dynamic.

```julia
ns = 10_000
S = generate_stacked_fourier_surrogates(Xres, ns)
lines!(axs[2][1], ttrend, S[1, :], label="surrogate residual")
lines!(axs[2][2], ttrend, S[ns+1, :], label="surrogate residual")
axislegend(axs[2][2], position = :lt)
fig
```

We then repeat the computation of indicators and their trends on the surrogate data, by using the same syntax as for the original residual. Finally, we can compute and visualize the percentile represented by the original time series with respect to the surrogate data. This is done by running:

```julia
surrogate_ar1 = slide_estimator( S, p_window_indctr, ar1_whitenoise, window )
surrogate_ar1_trend = slide( surrogate_ar1, tindctr, p_window_idtrend, ridge_regression_slope, window)
ar1_psignificance = percentile_significance(origin_ar1_trend, surrogate_ar1_trend, ns, nx)

for j in 1:nx
    lines!(axs[5][j], tidtrend, ar1_psignificance[j, :], label=L"$p$")
    hlines!(axs[5][j], [0.95], label=L"$95$th percentile", color = :red)
end
fig
```

A single transition indicator is often not robust enough to avoid false positives. In other words, the trend significance of the AR1 regression coefficient might lead to alarms when no transition is coming. Therefore, it is common to compute at least two transition indicators. In the next section we show how to compute arbitrarily many of them with a single line of code.

## [The fast-forward way](@id ffway)

In many cases, one does not need to perform each step of the computation. Therefore, a method is provided running all the computations given a matrix of residual time series:

```julia
sol = indicate_transition(
    ttrend,
    Xres,
    p_window_indctr,
    p_window_idtrend,
    ns,
    [ar1_whitenoise, var, lfps],
)
```

To plot the results, we run:

```julia

function plot_sol_doublewell(ttrend, Xtrend, sol, nx)

    nrows = 4
    fig = Figure( resolution = (1500, 1500) )
    axs = [[Axis(
        fig[i,j],
        xlabel = ( i==4 ? "Time" : " " ),
        xminorticks = 0:5:100,
        xminorgridvisible = true,
        ) for j in 1:2] for i in 1:4]

    for j in 1:nx

        lines!(axs[1][j], ttrend, Xtrend[j, :], label = L"trend $\,$")
        vlines!(
            axs[1][j],
            sol.predictor_time[j],
            label = L"warnings $\,$",
            color = :red,
            linestyle = :dash,
            linewidth = 1.0,
        )

        lines!(
            axs[2][j],
            sol.tindctr,
            sol.reference_ti["ar1_whitenoise"][j, :],
            label = L"AR1 $\,$",
        )
        lines!(
            axs[2][j],
            sol.tidtrend,
            sol.significance["ar1_whitenoise"][j, :], 
            label = L"AR1 significance $\,$",
        )

        lines!(
            axs[3][j],
            sol.tindctr,
            sol.reference_ti["var"][j, :] ./ maximum(sol.reference_ti["var"][j, :]),
            label = L"normalized variance $\,$",
        )
        lines!(
            axs[3][j],
            sol.tidtrend,
            sol.significance["var"][j, :], 
            label = L"normalized variance significance $\,$",
        )

        lines!(axs[4][j],
            sol.tindctr,
            sol.reference_ti["lfps"][j, :],
            label = L"LFPS $\,$",
        )
        lines!(
            axs[4][j],
            sol.tidtrend,
            sol.significance["lfps"][j, :], 
            label = L"LFPS significance $\,$",
        )

        [ylims!(axs[i][j], (-.1, 1.1)) for i in 2:4]
    end
    [axislegend(axs[i][2], position = :lt) for i in 1:4]

    return fig
end

fig_ffway = plot_sol_doublewell(ttrend, Xtrend, sol, nx)
```

Here the default settings were used to lump the percentile significance of the indicator trends (several continuous values) into a warning (a single binary value). Hence, a warning is only output when **all** the listed transition indicators display a percentile significance above 95%. If one wants a warning to be output when **at least two** indicators display a trend with a percentile significance above 99%, this can be done by running:

```julia
sol = indicate_transition(
    ttrend,
    Xres,
    p_window_indctr,
    p_window_idtrend,
    ns,
    [ar1_whitenoise, var, lfps],
    min_num_indicators = 3,
)
```

## GPU parallelisation

In the last decade, high-performance computing has made a drastic turn towards parallelism. In this context, Graphics Processing Units (GPUs) have become a popular solution as they offer a good performance vs. cost and benefit from increasingly user-friendly interfaces.

The present package allows to run the most expensive part of the computation on a GPU. To do so, one simply needs to convert the array containing the residuals into a GPU array:

```julia
Xres_gpu = CuArray(Xres)
sol = indicate_transition(
    ttrend,
    Xres_gpu,
    p_window_indctr,
    p_window_idtrend,
    ns,
    [ar1_whitenoise, var, lfps],
)
```

The present section only aims to show the addition of a single line of code to run computation on the GPU. For a thorough analysis of the speed increase, refer to [Benchmarks](@ref benchmarks). 

!!! note
For now, the present package only supports NVIDIA GPUs as it merely relies on [CUDA.jl](https://cuda.juliagpu.org/stable/).
!!!