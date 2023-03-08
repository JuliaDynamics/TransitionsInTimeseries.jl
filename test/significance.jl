using TimeseriesSurrogates, TransitionIndicators, Test, Random, Statistics

curry(f, y) = x -> f(x, y)

function generate_results()
    t = collect(0.0:100_000.0)
    θ = rand()
    x = AR1(length(t), rand(), θ, Random.default_rng())

    p = SignificanceHyperParams(
        n_surrogates = 1_000,
        wv_indicator_width = 10,
        wv_evolution_width = 10,
    )
    evolution_metric = precomputed_ridge_slope(p)
    res = analyze_indicators(t, x, var, evolution_metric, p)
    return res
end

function count_percentile(res, sig)
    significance = vec(measure_significance(res, sig))
    top5 = sum(significance .> 1.0)
    bottom5 = sum(significance .< -1.0)
    return bottom5 + top5
end

@testset "measuring significance w.r.t. to n surrogates" begin
    # Generate long time series to check if statistics are well-behaved
    res = generate_results()

    # We look at 95-th percentile (either symmetric or only one-sided).
    # The significance of our AR1 process vs. AR1 surrogates should therefore be
    # positive in [4, 6]% of the cases.
    tol = intround.(length(res.X_evolution[1, :, 1]) .* [0.04, 0.06])

    # Test if confidence_interval() gives significance within tolerance
    @test count_percentile(res, confidence_interval) in tol[1]:tol[2]

    # Test if symmetric_nqd() gives significance within tolerance
    symmetric_nqd(x, s) = normalized_percentile(x, s, symmetric = true)
    @test count_percentile(res, symmetric_nqd) in tol[1]:tol[2]

    # Test if asymmetric_nqd() gives significance within tolerance
    asymmetric_nqd(x, s) = normalized_percentile(x, s)
    @test count_percentile(res, asymmetric_nqd) in tol[1]:tol[2]

    # Test if which_percentile() gives significance within tolerance
    # which_percentile() is not normed and therefore needs special treatment
    significance = measure_significance(res, which_percentile)
    significance = vec(significance)
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

    nid5 = normalized_percentile(3.0, s, p = 0.95, symmetric = true)
    nid50 = normalized_percentile(50.5, s, p = 0.95, symmetric = true)
    nid95 = normalized_percentile(98.0, s, p = 0.95, symmetric = true)

    @test isapprox(nid5, -1)
    @test isapprox(nid50, 0)
    @test isapprox(nid95, 1)
end