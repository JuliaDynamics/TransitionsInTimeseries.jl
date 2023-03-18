using TransitionIndicators, Test, Random, Statistics

function generate_results()
    t = collect(0:10_000)
    θ = rand()
    x = AR1(length(t), rand(), θ, Random.default_rng())

    indicators = [var, ar1_whitenoise]
    ind_conf = IndicatorsConfig(indicators; width = 100, stride = 1)

    change_width = 10
    change_metric = RidgeRegression(t, change_width)
    sig_conf = SignificanceConfig(change_metric; width = change_width, n_surrogates = 5_000)

    return indicators_analysis(x, ind_conf, sig_conf)
end

@testset "measuring significance w.r.t. to n surrogates" begin
    # Generate long time series to check if statistics are well-behaved
    res = generate_results()

    # We look at 95-th quantile. The significance of AR1 process vs. AR1 surrogates
    # should therefore be positive in about 5% of the cases.
    # Here we set the acceptable range to be [4, 6]%.
    tol = size(res.x_change, 1) .* [0.04, 0.06]
    up95quantile = Quantile(0.95, :up)
    sig = indicators_significance(res, up95quantile)
    @test tol[1] < sum(sig[:,1]) < tol[2]
    @test tol[1] < sum(sig[:,2]) < tol[2]

    # We look at (upper) 2σ confidence intervall. The significance of AR1 process
    # vs. AR1 surrogates should therefore be positive in about 2.5% of the cases.
    # Here we set the acceptable range to be [2, 3]%.
    tol = size(res.x_change, 1) .* [0.02, 0.03]
    up95confintervall = Sigma(2, :up)
    sig = indicators_significance(res, up95confintervall)
    @test tol[1] < sum(sig[:,1]) < tol[2]
    @test tol[1] < sum(sig[:,2]) < tol[2]
end

@testset "significant quantile" begin
    s = 1:100

    updown95quantile = Quantile(0.95, :updown)
    down95quantile = Quantile(0.95, :down)
    up95quantile = Quantile(0.95, :up)

    @test significant(5, s, updown95quantile)
    @test !significant(6, s, updown95quantile)
    @test !significant(95, s, updown95quantile)
    @test significant(96, s, updown95quantile)

    @test significant(5, s, down95quantile)
    @test !significant(6, s, down95quantile)
    @test !significant(95, s, down95quantile)
    @test !significant(96, s, down95quantile)

    @test !significant(5, s, up95quantile)
    @test !significant(6, s, up95quantile)
    @test !significant(95, s, up95quantile)
    @test significant(96, s, up95quantile)
end

@testset "significant confidence intervall" begin
    s = randn(10_000)
    updown95confintervall = Sigma(2, :updown)
    down95confintervall = Sigma(2, :down)
    up95confintervall = Sigma(2, :up)

    @test significant(-2.1, s, updown95confintervall)
    @test !significant(-1.9, s, updown95confintervall)
    @test !significant(1.9, s, updown95confintervall)
    @test significant(2.1, s, updown95confintervall)

    @test significant(-2.1, s, down95confintervall)
    @test !significant(-1.9, s, down95confintervall)
    @test !significant(1.9, s, down95confintervall)
    @test !significant(2.1, s, down95confintervall)

    @test !significant(-2.1, s, up95confintervall)
    @test !significant(-1.9, s, up95confintervall)
    @test !significant(1.9, s, up95confintervall)
    @test significant(2.1, s, up95confintervall)
end