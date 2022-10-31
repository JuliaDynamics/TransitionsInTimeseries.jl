# Example: double-well system

## Introduction

On of the simplest system displaying bi-stability is the double-well potential. By taking the negative gradient of such 1D potential, one obtains following normal form:

```math
\dot{x} = -x^3 + x + F(t) + n(t)
```

with $\dot{x}$ the time derivative of the state $x$, $F$ the forcing and $n$ a source of additive white noise with variance $\sigma^2$. For the autonomous case, the system displays stable equilibrium points in $x=-1$, $x=1$ and an unstable one in $x=0$. By setting $x(t=0) = -1$ and applying a slow linear drift in $F$, one can artificially generate a transition, which is expected to be detectable thanks to TIs.

A common problem of TIs is to deliver false positives, i.e. prediction of a transition ahead even when none is about to happen. To investigate this aspect in parallel, we study a linear model also displaying an equilibrium point at $x=-1$ for the autonomous case and defined as:

```math
\dot{x} = \lambda (x+1) + F(t) + n(t).
```

As linear models do not allow critical transitions to occur, they provide a baseline to check for false positives. Putting these equations into code results in:

```julia
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

Let's now define the models with initial condition in their common equilibrium point $x=-1$ and generate forced time series.

```julia
using DifferentialEquations

T = Float32
dt = 5f-2
t = 0f0:dt:100f0
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

nrows, ncols = 5, nx
fig = Figure( resolution = (1500, 1500) )
axs = [[Axis(
    fig[i,j],
    title = ( i==1 ? labels[j] : " "),
    xlabel = ( i==nrows ? "Time" : " "),
    xminorticks = tspan[1]:5:tspan[2],
    xminorgridvisible = true,
    ) for j in 1:ncols] for i in 1:nrows]

for j in 1:nx
    lines!(axs[1][j], t, X[j, :], label = L"data $\,$")
    ylims!(axs[1][j], (-1.2, 1.5))
fig
```

## Smoothing and de-trending

Critical slowing down is possibly observed in the response of  the system to noise. To study this without taking the slower, deterministic dynamic into account, the time series is filtered. Here some of the many options are:

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

An extensive list of kernels can be found [here](https://en.wikipedia.org/wiki/Kernel_(statistics)). They are implemented in the present package and can be found in the [`API reference`](@ref api_ref). In the present case, a centered window is applied on the data $x$, but the user can choose whether they prefer to use e.g. a left-window. Furthermore, the half-width $N$ of the window as well as the stride with which the operation should be applied are to be determined by the user. An exemplary code is here provided:

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
    lines!(axs[2][i], ttrend, Xres[i, :], label = L"residual $\,$")
end

axislegend(axs[1][2], position = :lt)
fig
```
## The step-by-step way

From here on, two paths open for the user: either do the computation step by step and keep a large flexibility to process data in between, or go the fast-forward way to directly obtain the transition indicators, their trends, their significance and associated potential warning for transition.

Whatever way you prefer, we recommend to read the step-by-step way to understand what happens under the hood of [The fast-forward way](@ref ffway).

### Estimating the transition indicators

Now that the residuals are computed, one can estimate transition indicators. For this we first need to define windowing parameters and compute the associated time vector:

```julia
# Windowing parameters for TI estimation.
T_indctr_wndw = T(5e-1)
T_indctr_strd = T(dt)
p_window_indctr = get_windowing_params([dt, T_indctr_wndw, T_indctr_strd])

window = left_wndw
tindctr = trim_wndw(ttrend, p_window_indctr, window)
```

The implemented functions to estimate transition indicators is provided in the [API reference](@ref api_indicators). By defining a list with such functions, one can easily specify which transition indicators are to be computed. The computation itself is performed via the [`slide_estimator`](@ref) function, which slides the windowed computation over the time series.

```julia
TIlist = [var, ar1_whitenoise, lfps]
nti = length(TIlist)

# pre-allocate an array for TI computation.
TIref = zeros(T, nx, length(tindctr), nti)

# labels associated with each TI for plotting.
TIlabels = [
    L"variance $\,$",
    L"AR1 coefficient $\,$",
    L"LFPS $\,$"]

for i in 1:nti
    TIref[:, :, i] = slide_estimator(Xres, p_window_indctr, TIlist[i], window)
    
    for j in 1:nx
        lines!(axs[2+i][j], tindctr, TIref[j, :, i], label = TIlabels[i])
    end
end
[ylims!(axs[i][j], (-.1, 1.1)) for i in 3:nrows, j in 1:ncols]

fig
```

### Estimating indicator trends

To indicate a transition, one is primarly interested in the trend of an indicator rather than in its absolute value. Here again, the computation requires the definition of a certain window over which the trend computation is performed:

```julia
T_idtrend_wndw = T(2e0)
T_idtrend_strd = T(1e0)
dt_id = round(tindctr[2] - tindctr[1]; digits = 4)

p_window_idtrend = get_windowing_params([dt_id, T_idtrend_wndw, T_idtrend_strd])
tidtrend = trim_wndw( tindctr, p_window_idtrend, window )
```

Once this is done, one can easily compute the trends of the transition indicators by using one of the methods listed in the [API reference](@ref api_trends)

```julia
TItrend_ref = zeros(T, nx, length(tidtrend), nti)

for i in 1:nti
    for j in 1:nx
        TItrend_ref = slide_idtrend
```

### Significance of indicator trends


## [The fast-forward way](@id ffway)

Now that the residuals are computed, one can obtain transition indicators, their trends and their significance by running:

```julia
# Number of surrogates per time series.
ns = 10_000



# Windowing params for TI trend estimation.
T_idtrend_wndw = T(2e0)
T_idtrend_strd = T(1e0)
p_window_idtrend = get_windowing_params([dt, T_idtrend_wndw, T_idtrend_strd])

sol = identify_transition(
    ttrend,
    Xres,
    p_window_indctr,
    p_window_idtrend,
    ns,
    [ar1_whitenoise, var, lfps],
)
```

