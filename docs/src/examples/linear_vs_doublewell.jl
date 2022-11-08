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

######################

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

######################

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

######################

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

######################

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

######################

T_idtrend_wndw = T(2e0)
T_idtrend_strd = T(1e0)
dt_id = round(tindctr[2] - tindctr[1]; digits = 4)
p_window_idtrend = get_windowing_params([dt_id, T_idtrend_wndw, T_idtrend_strd])

tidtrend = trim_wndw( tindctr, p_window_idtrend, window )
origin_ar1_trend = slide_idtrend( origin_ar1, tindctr, p_window_idtrend, ridge_regression_slope, window)

for j in 1:nx
    lines!(axs[4][j], tidtrend, origin_ar1_trend[j, :], label = L"trend of AR1 coefficient $\,$")
end
fig

######################

ns = 10_000
S = generate_stacked_fourier_surrogates(Xres, ns)
lines!(axs[2][1], ttrend, S[1, :], label="surrogate residual")
lines!(axs[2][2], ttrend, S[ns+1, :], label="surrogate residual")
axislegend(axs[2][2], position = :lt)
fig

######################

surrogate_ar1 = slide_estimator( S, p_window_indctr, ar1_whitenoise, window )
surrogate_ar1_trend = slide_idtrend( surrogate_ar1, tindctr, p_window_idtrend, ridge_regression_slope, window)
ar1_psignificance = percentile_significance(origin_ar1_trend, surrogate_ar1_trend, ns, nx)

for j in 1:nx
    lines!(axs[5][j], tidtrend, ar1_psignificance[j, :], label=L"$p$")
    hlines!(axs[5][j], [0.95], label=L"$95$th percentile", color = :red)
end
fig

######################

sol = indicate_transition(
    ttrend,
    Xres,
    p_window_indctr,
    p_window_idtrend,
    ns,
    [ar1_whitenoise, var, lfps],
)

######################

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

######################


######################
