using TransitionIndicators
using CairoMakie

t, x_linear, x_nlinear = generate_test_data()
fig = Figure()
axs = [Axis(fig[i, j], xticklabelsvisible = i == 5 ? true : false) for i in 1:5, j in 1:2]
lines!(axs[1, 1], t, x_linear)
lines!(axs[1, 2], t, x_nlinear)
fig

#################################################

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

#################################################

using Statistics
p = init_metaanalysis_params()

indicators = [var, ar1_whitenoise]
evolution_metrics = precomputed_ridge_slope(p)

results_lin = analyze_indicators(t, lin_fluctuations, indicators, evolution_metrics, p)
results_nonlin = analyze_indicators(t, nlin_fluctuations, indicators, evolution_metrics, p)

results = [results_lin, results_nonlin]
indicator_number = 1
for j in eachindex(results)
    lines!(axs[3, j], results[j].t_indicator, results[j].X_indicator[1, :, indicator_number])
    lines!(axs[4, j], results[j].t_evolution, results[j].X_evolution[1, :, indicator_number])
end
fig

sig_lin = measure_significances(results_lin, normalized_confidence_intervall)
sig_nonlin = measure_significances(results_nonlin, normalized_confidence_intervall)
significances = [sig_lin, sig_nonlin]
for i in eachindex(indicators), j in eachindex(results)
    lines!(axs[5, j], results[j].t_evolution, significances[j][1, :, i])
end
fig








# @btime results_lin = analyze_indicators(t, lin_fluctuations, indicators, evolution_metrics, p)
# @btime results_nonlin = analyze_indicators(t, nlin_fluctuations, indicators, evolution_metrics, p)



# t, x_lin, x_nonlin = generate_test_data()
# fig, axs = init_plot()
# lines!(axs[1][1], t, x_lin, label = L"raw $\,$")
# lines!(axs[1][2], t, x_nonlin, label = L"raw $\,$")

# x_lin_smooth = running_mean(x_lin)
# x_nonlin_smooth = running_mean(x_nonlin)
# lines!(axs[1][1], t, x_lin_smooth, label = L"smooth $\,$")
# lines!(axs[1][2], t, x_nonlin_smooth, label = L"smooth $\,$")
# axislegend(axs[1][2], position = :rb)

# residual_linear = x_lin - x_lin_smooth
# residual_nonlinear = x_nonlin - x_nonlin_smooth
# lines!(axs[2][1], t, residual_linear)
# lines!(axs[2][2], t, residual_nonlinear)

# ###################################################

# function lfps(
#     x::AbstractVector{T};
#     q_lofreq::AbstractFloat=0.1,
# ) where {T<:Real}

#     power = abs.(rfft(x)) .^ 2
#     normed_power = power ./ sum(power)
#     lofreq_idx = Int(round(q_lofreq * length(normed_power)))
#     return sum( normed_power[1:lofreq_idx] )
# end

# p = init_metaanalysis_params(
#     n_surrogates = 10_000,
#     surrogate_method = RandomFourier(),
#     wv_indicator_width = 50,
#     wv_indicator_stride = 5,
#     wv_evolution_width = 40,
#     wv_evolution_stride = 5,
# )
# m = precompute_ridge_slope(1.0:p.wv_evolution_width+1, lambda = 0.1)
# indicators = [var, ar1_whitenoise, lfps]
# evolution_metrics = curry(precomputed_ridge_slope, m)
# # evolution_metrics = [kendalltau, curry(precomputed_ridge_slope, m), spearman]
# ni = length(indicators)

# results_lin = analyze_indicators(
#     t, residual_linear, indicators, evolution_metrics, p )
# results_nonlin = analyze_indicators(
#     t, residual_nonlinear, indicators, evolution_metrics, p )
# results = [results_lin, results_nonlin]
# colors = [:blue, :orange]
# for i in 1:2, m in 1:2
#     lines!(axs[3][m], results[m].t_indicator, results[m].X_indicator[1, :, i])
#     lines!(axs[4][m], results[m].t_evolution, results[m].X_evolution[1, :, i])
# end


# sig_lin = measure_significances(results_lin, normalized_confidence_intervall)
# sig_nonlin = measure_significances(results_nonlin, normalized_confidence_intervall)
# significances = [sig_lin, sig_nonlin]
# for i in 1:ni, m in 1:2
#     lines!(axs[5][m], results[m].t_evolution, significances[m][1, :, i])
# end

# fig