using TransitionsInTimeseries, Test, Random, Distributions

@testset "ridge regression" begin
    d = Normal()                            # define normal distribution
    m, p = 10 .* rand(Xoshiro(124), d, 2)   # slope and offset parameters
    t = 0:0.1:100
    x = m .* t .+ p
    x_noisy = x + 0.1 .* rand(d, length(x))

    rr = precompute(RidgeRegressionSlope(), t)
    m_hat = rr(x)
    m_hat_noisy = rr(x_noisy)
    @test isapprox(m_hat, m, atol = 1e-5)
    @test isapprox(m_hat_noisy, m, atol = 1e-2)
end

@test "correlations" begin
    x = 0:0.01:2π
    y = sin.(x)
    @test spearman(x) == 1
    # TODO: expand this when kendalltau fixed
end