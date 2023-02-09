using TransitionIndicators, Test, Statistics

@testset "computing normalized quantile distance" begin
    s = 1.0:100.0
    S = repeat(s, outer = (1, 20))
    symmetric_quantile_idx = quantile_idx(s, q = 0.95, symmetric = true)
    asymmetric_quantile_idx = quantile_idx(s, q = 0.95, symmetric = false)
    @test symmetric_quantile_idx == (5, 95)
    @test asymmetric_quantile_idx == (1, 95)

    nid5 = normalized_quantile_distance(5.0, s, q=0.95, symmetric = true)
    nid50 = normalized_quantile_distance(50.0, s, q=0.95, symmetric = true)
    nid95 = normalized_quantile_distance(95.0, s, q=0.95, symmetric = true)

    @test isapprox(nid5, -1)
    @test isapprox(nid50, 0)
    @test isapprox(nid95, 1)
end

@testset "measuring significance w.r.t. to n surrogates" begin
    # Generate long time series to check if statistics are well-behaved
    t = collect(0.0:10_000.0)
    x = AR1(length(t), 0.0, -0.5)
    
    wv_indicator_width = 10
    wv_indicator_stride = 2
    wv_evolution_width = 10
    wv_evolution_stride = 2

    M = precompute_ridge_slope(1:wv_evolution_width+1)
    curry(f, y) = x -> f(x, y)
    trend_results = indicator_evolution(
        t,
        x,
        Statistics.var,
        1_000,
        RandomFourier(),
        curry(precomputed_ridge_slope, m),
        wv_indicator_width,
        wv_indicator_stride,
        wv_evolution_width,
        wv_evolution_stride,
    )

    # The significance of our AR1 process vs. AR1 surrogates should be
    # positive in [4, 6]% of the cases.
    tol = intround.(length(trend_results.x_evolution) .* [0.04, 0.06])
    tol_less = 0.01 * length(trend_results.x_evolution)

    # Test gaussian_quantile()
    significance = measure_significance(trend_results, gaussian_quantile)
    gq_top5 = sum(significance .> 1.0)
    gq_bottom5 = sum(significance .< -1.0)
    @test (gq_top5 + gq_bottom5) in tol[1]:tol[2]

    # Test symmetric normalized_quantile_distance()
    symmetric_nqd(x, s) = normalized_quantile_distance(x, s, symmetric = true)
    significance = measure_significance(trend_results, symmetric_nqd)
    snqd_top5 = sum(significance .> 1.0)
    snqd_bottom5 = sum(significance .< -1.0)
    @test snqd_top5 in tol[1]:tol[2]
    @test snqd_bottom5 in tol[1]:tol[2]

    # Test asymmetric normalized_quantile_distance()
    asymmetric_nqd(x, s) = normalized_quantile_distance(x, s)
    significance = measure_significance(trend_results, asymmetric_nqd)
    anqd_top5 = sum(significance .> 1.0)
    anqd_bottom5 = sum(significance .< -1.0)
    @test anqd_top5 in tol[1]:tol[2]
    @test anqd_bottom5 < tol_less

    # Test which_quantile()
    significance = measure_significance(trend_results, which_quantile)
    top5 = sum(significance .> 0.95)
    bottom5 = sum(significance .< 0.05)
    @test top5 in tol[1]:tol[2]
    @test bottom5 in tol[1]:tol[2]
end