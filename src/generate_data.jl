"""

    generate_test_data()

Generate prototypical data from a linear and a double-fold model to test some indicators.
Keyword arguments:
- `T::Type`, the output type
- `dt::Real` = 5e-2, the fixed time step for the output
- `t_end::Real` = 50.0, the final integration time
- `x0::Real` = -1.0, the initial condition
- `alpha::Real` = 0.02, the slope of forcing ramp
- `lambda::Real` = -1.0, the resotring rate of linear system
- `sigma::Real` = 0.02, the standard-deviation of the white noise
"""
function generate_test_data(;
    T::Type = Float64,
    dt::Real = 5e-2,
    t_end::Real = 50.0,
    x0::Real = -1.0,
    alpha::Real = 0.02,
    lambda::Real = -1.0,
    sigma::Real = 0.02,
)

    t = collect(T(0):T(dt):T(t_end))
    tspan = extrema(t)
    x0 = [T(x0)]
    pmodel = Dict(
        "α" => alpha,
        "λ" => lambda,
        "σ" => sigma,
    )
    models = [f_linear, f_doublewell]
    nx = length(models)
    X = zeros(T, nx, length(t))

    for (f, i) in zip(models, 1:nx)
        prob = SDEProblem(f, g_whitenoise, x0, tspan, pmodel)
        sol = solve(prob, EM(), dt=dt)
        X[i, :] .= vcat(sol.u...)[1:end]
    end
    return t, X[1, :], X[2, :]
end

# Deterministic part of linear system.
function f_linear(dx, x, p, t)
    dx[1] = p["λ"] * (x[1] + 1) + ramp_forcing(p, t)
end

# Deterministic part of double-well system.
function f_doublewell(dx, x, p, t)
    dx[1] = (- x[1]^3 + x[1]) + ramp_forcing(p, t)
end

# Stochastic part of both system for noise.
function g_whitenoise(dx, x, p, t)
    dx[1] = p["σ"]
end

# Ramp forcing
function ramp_forcing(p, t)
    return p["α"] * t
end

