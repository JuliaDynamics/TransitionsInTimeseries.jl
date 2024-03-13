using TransitionsInTimeseries, Test, Random, Distributions

@testset "ridge regression" begin
    d = Normal()                            # define normal distribution
    m, p = 10 .* rand(Xoshiro(124), d, 2)   # slope and offset parameters
    t = 0:1:100
    x = m .* t .+ p
    x_noisy = x + 0.1 .* rand(d, length(x))

    rr = precompute(RidgeRegressionSlope(), t)
    m_hat = rr(x)
    m_hat_noisy = rr(x_noisy)
    @test isapprox(m_hat, m, atol = 1e-5)
    @test isapprox(m_hat_noisy, m, atol = 1e-2)

    rr = RidgeRegressionSlope()
    m_hat = rr(x)
    m_hat_noisy = rr(x_noisy)
    @test isapprox(m_hat, m, atol = 1e-5)
    @test isapprox(m_hat_noisy, m, atol = 1e-2)
end

@testset "correlations" begin
    x = collect(1:5)
    y = [1, 2, 3, 5, 4]
    @test spearman(x) == 1
    @test kendalltau(x) == 1
    @test isapprox(spearman(y), cov(x, y) / (std(x) * std(y)))
    @test kendalltau(y) == 0.8
end