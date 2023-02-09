using Random, Distributions, TransitionIndicators, Test

@testset "checking slope estimation" begin
    d = Normal()                # define normal distribution
    m, p = 10 .* rand(d, 2)     # slope and offset parameters
    x = 0:0.1:10
    y = m .* x .+ p             # define affine function
    y_noisy = y + 0.1 .* rand(d, length(x))
    m_est, p_est = ridge(x, y)
    m_est_noisy, p_est_noisy = ridge(x, y_noisy)
    @test isapprox(m, m_est)
    @test isapprox(p, p_est)
    @test isapprox(m, m_est_noisy, atol = 1e-1)
    @test isapprox(p, p_est_noisy, atol = 1e-1)
end

@testset "sliding trend estimation over indicator" begin
    n = 101
    t = collect(1.0:n)
    x = copy(t)
    
    indicator_wv_width = 10
    indicator_wv_stride = 2
    metric_wv_width = 10
    metric_wv_stride = 2
    slope_ground_truth = fill(indicator_wv_stride, length(slope_ts))

    trend_results = indicator_evolution(
        t,
        x,
        mean,
        100,
        RandomFourier(),
        ridge_slope,
        indicator_wv_width,
        indicator_wv_stride,
        metric_wv_width,
        metric_wv_stride,
    )
    @test isapprox(trend_results.x_evolution, slope_ground_truth)

end