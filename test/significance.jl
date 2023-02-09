using TransitionIndicators, Test, Statistics

@testset "measuring significance w.r.t. to n surrogates" begin
    # Generate long time series to check if statistics are well-behaved
    t = collect(0.0:10_000.0)
    x = AR1(length(t), 0.0, -0.5)
    
    wv_indicator_width = 10
    wv_indicator_stride = 2
    wv_evolution_width = 10
    wv_evolution_stride = 2

    m = precompute_ridge_slope(1:wv_evolution_width+1)
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

    # We look at 95-th percentile (either symmetric or only one-sided).
    # The significance of our AR1 process vs. AR1 surrogates should therefore be
    # positive in [4, 6]% of the cases.
    tol = intround.(length(trend_results.x_evolution) .* [0.04, 0.06])

    # Test gaussian_percentile()
    significance = measure_significance(trend_results, gaussian_percentile)
    gq_top5 = sum(significance .> 1.0)
    gq_bottom5 = sum(significance .< -1.0)
    @test (gq_top5 + gq_bottom5) in tol[1]:tol[2]

    # Test symmetric normalized_percentile_distance()
    symmetric_nqd(x, s) = normalized_percentile_distance(x, s, symmetric = true)
    significance = measure_significance(trend_results, symmetric_nqd)
    snqd_top5 = sum(significance .> 1.0)
    snqd_bottom5 = sum(significance .< -1.0)
    @test (snqd_top5 + snqd_bottom5) in tol[1]:tol[2]

    # Test asymmetric normalized_percentile_distance()
    asymmetric_nqd(x, s) = normalized_percentile_distance(x, s)
    significance = measure_significance(trend_results, asymmetric_nqd)
    anqd_top5 = sum(significance .> 1.0)
    anqd_bottom5 = sum(significance .< -1.0)
    @test (anqd_top5 + anqd_bottom5) in tol[1]:tol[2]

    # Test which_percentile()
    significance = measure_significance(trend_results, which_percentile)
    top5 = sum(significance .> 0.975)
    bottom5 = sum(significance .< 0.025)
    @test (top5 + bottom5) in tol[1]:tol[2]
end

@testset "computing normalized percentile distance" begin
    s = 1.0:100.0
    S = repeat(s, outer = (1, 20))
    symmetric_percentile_idx = percentile_idx(s, p = 0.95, symmetric = true)
    asymmetric_percentile_idx = percentile_idx(s, p = 0.95, symmetric = false)
    @test symmetric_percentile_idx == (3, 98)
    @test asymmetric_percentile_idx == (1, 95)

    nid5 = normalized_percentile_distance(3.0, s, p = 0.95, symmetric = true)
    nid50 = normalized_percentile_distance(50.5, s, p = 0.95, symmetric = true)
    nid95 = normalized_percentile_distance(98.0, s, p = 0.95, symmetric = true)

    @test isapprox(nid5, -1)
    @test isapprox(nid50, 0)
    @test isapprox(nid95, 1)
end