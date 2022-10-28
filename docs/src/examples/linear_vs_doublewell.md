# Example double-well

## Introduction

On of the simplest system displaying bi-stability is the double-well potential. By taking the negative gradient of such 1D potential, one obtains following normal form:

$$
\dot{x} = -x^3 + x + F(t) + n(t)
$$

with $\dot{x}$ the time derivative of the state $x$, $F$ the forcing and $n$ a source of additive white noise with variance $\sigma^2$. For the autonomous case, the system displays stable equilibrium points in $x=-1$, $x=1$ and an unstable one in $x=0$. By setting $x(t=0) = -1$ and applying a slow linear drift in $F$, one can artificially generate a transition, which is expected to be detectable thanks to TIs.

A common problem of TIs is to deliver false positives, i.e. prediction of a transition ahead even when none is about to happen. To investigate this aspect in parallel, we study a linear model also displaying an equilibrium point at $x=-1$ for the autonomous case and defined as:

$$
\dot{x} = \lambda (x+1) + F(t) + n(t).
$$

As linear models do not allow critical transitions to occur, they provide a baseline to check for false positives. 

## Generating data of critical transition

Let's now define the models with initial condition in their common equilibrium point $x=-1$ and generate forced time series.

```julia
using TransitionIdentifiers
using DifferentialEquations, CairoMakie

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

T = Float32
dt = 5f-2
t = 0f0:dt:100f0
tspan = extrema(t)
x0 = [-1f0]

pmodel = Dict(
    "σ" => T(0.02),
    "λ" => T(-1),
    "α" => T(0.01),
)

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

$$
X = \begin{pmatrix}
& x^{(1)}_1 & x^{(1)}_2 & \dots & x^{(1)}_{n_t- 1} & x^{(1)}_{n_t} & \\
& \\
& x^{(2)}_1 & x^{(2)}_2 & \dots & x^{(2)}_{n_t- 1} & x^{(2)}_{n_t} &
\end{pmatrix},
$$

where $x^{(j)}_k$ the state of model $j$ (1=linear, 2=double-well) at time step $k$. Let us briefly visualise the time series:

```julia
nrows, ncols = 5, nx
fig = Figure( resolution = (1500, 1500) )
axs = [[Axis(
    fig[i,j],
    title = ( i==1 ? labels[j] : " "),
    xlabel = ( i==nrows ? "Time" : " "),
    xminorticks = tspan[1]:5:tspan[2],
    xminorgridvisible = true,
    ) for j in 1:ncols] for i in 1:nrows]

[lines!(axs[1][i], t, X[i, :], label = L"data $\,$") for i in 1:ncols]
[ylims!(axs[1][i], (-1.2, 1.5)) for i in 1:ncols]
fig
```

## Smoothing and de-trending

## Estimating indicators



## Estimating indicator trends

## Significance of indicator trends
