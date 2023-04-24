using TransitionIndicators, Test, Random, TimeseriesSurrogates
using Distributions

# @testset "checking slope estimation" begin
#     d = Normal()                # define normal distribution
#     m, p = 10 .* rand(d, 2)     # slope and offset parameters
#     x = 0:0.1:10
#     y = m .* x .+ p             # define affine function
#     y_noisy = y + 0.1 .* rand(d, length(x))
#     rr = RidgeRegression(x)
#     m_est = rr(y)
#     m_est_noisy = rr(y_noisy)
#     @test isapprox(m, m_est)
#     @test isapprox(m, m_est_noisy, atol = 1e-1)
# end

d = Normal()                # define normal distribution
m, p = 10 .* rand(d, 2)     # slope and offset parameters
t = 0:0.1:100
x = m .* t .+ p             # define affine function
noise = 0.1 .* rand(d, length(x))
x_noisy = x + noise

iparams = IndicatorsParams(q_lofreq = 0.1)
indconfig = IndicatorsConfig(t, [ar1_whitenoise, LowfreqPowerSpectrum],
    width = 100, stride = 2)

cmparams = ChangeMetricsParams(lambda = 0.0)
sigconfig = SignificanceConfig(
    indconfig,
    [RidgeRegression, kendalltau],
    change_metrics_params = cmparams,
    width = 10,
    stride = 1,
)

lfps_test = indconfig.indicators[2](noise[1:indconfig.width])
x_indicator = x[indconfig.width:indconfig.stride:end]
sigconfig.change_metrics[1](x_indicator[1:sigconfig.width])

# @testset "sliding trend estimation over indicator" begin
#     n = 1001
#     t = collect(1.0:n)
#     x = copy(t)
#     p = SignificanceHyperParams(n_surrogates = 100)
#     rr = RidgeRegression(t, p.wv_evolution_width)
#     res = analyze_indicators(t, x, [mean, var], rr, p)

#     # The trend of mean(windowview) is the stride for x=t
#     meantrend_ground_truth = fill(p.wv_indicator_stride, length(res.t_evolution))
#     # The trend of var(windowview) is 0 for x any affine function of t.
#     vartrend_ground_truth = fill(0.0, length(res.t_evolution))
#     @test isapprox(res.X_evolution[1, :, 1], meantrend_ground_truth)
#     @test isapprox(res.X_evolution[1, :, 2], vartrend_ground_truth, atol = 1e-12)
# end
