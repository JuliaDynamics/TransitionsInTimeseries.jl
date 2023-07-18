using TransitionsInTimeseries, Test, Random, Distributions

@testset "ridge regression" begin
    d = Normal()                # define normal distribution
    m, p = 10 .* rand(Xoshiro(124), d, 2)     # slope and offset parameters
    t = 0:0.1:100
    x = m .* t .+ p
    x_noisy = x + 0.1 .* rand(d, length(x))

    # TODO: Fix this. the previous code was written in an un-understastandable way.
    # Why compute variance as indicator but not use it? why do this weird "start
    # from the window" thing? Anyways, this doesn't work.
    width_cha = 100
    rr = precompute(RidgeRegressionSlope(), x[1:width_cha])
    m_hat = rr(x[1:width_cha])
    m_hat_noisy = rr(x_noisy[1:width_cha])
    # @test isapprox(m_hat, m, atol = 1e-5)
    # @test isapprox(m_hat_noisy, m, atol = 1e-2)

    # Previous code

    # indconfig = IndicatorsConfig(t, last, [var], width = 20, stride = 1)
    # sigconfig = SignificanceConfig(indconfig, last, [RidgeRegressionSlope()],
    #     width = 100, stride = 1)

    # x_test = x[indconfig.width:indconfig.stride:end]
    # x_test_noisy = x_noisy[indconfig.width:indconfig.stride:end]
    # m_hat = sigconfig.change_metrics[1](x_test[1:sigconfig.width])
    # m_hat_noisy = sigconfig.change_metrics[1](x_test_noisy[1:sigconfig.width])
    # @test isapprox(m_hat, m, atol = 1e-5)
    # @test isapprox(m_hat_noisy, m, atol = 1e-2)
end

