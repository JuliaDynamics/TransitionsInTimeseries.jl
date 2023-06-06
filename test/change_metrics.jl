using TransitionsInTimeseries, Test, Random, Distributions

@testset "ridge regression" begin
    d = Normal()                # define normal distribution
    m, p = 10 .* rand(d, 2)     # slope and offset parameters
    t = 0:0.1:100
    x = m .* t .+ p             # define affine function
    x_noisy = x + 0.1 .* rand(d, length(x))

    indconfig = IndicatorsConfig(t, last, [var], width = 20, stride = 1)
    sigconfig = SignificanceConfig(indconfig, last, [RidgeRegressionSlope()],
        width = 100, stride = 1)

    x_test = x[indconfig.width:indconfig.stride:end]
    x_test_noisy = x_noisy[indconfig.width:indconfig.stride:end]
    m_hat = sigconfig.change_metrics[1](x_test[1:sigconfig.width])
    m_hat_noisy = sigconfig.change_metrics[1](x_test_noisy[1:sigconfig.width])
    @test isapprox(m_hat, m, atol = 1e-5)
    @test isapprox(m_hat_noisy, m, atol = 1e-2)
end


@testset "sliding ridge regression" begin
    n = 1001
    t = collect(1.0:n)
    x = copy(t)

    indconfig = IndicatorsConfig(t, last, [mean, var], width = 100, stride = 1)
    sigconfig = SignificanceConfig(indconfig, last,
        [RidgeRegressionSlope()], width = 100, stride = 1)
    res = indicators_analysis(t, x, indconfig, sigconfig)

    # The trend of mean(windowview) is the stride for x=t
    meantrend_ground_truth = fill(indconfig.stride, length(res.t_change))
    # The trend of var(windowview) is 0 for x any affine function of t.
    vartrend_ground_truth = fill(0.0, length(res.t_change))
    @test isapprox(res.x_change[:, 1], meantrend_ground_truth, atol = 1e-9)
    @test isapprox(res.x_change[:, 2], vartrend_ground_truth, atol = 1e-9)
end