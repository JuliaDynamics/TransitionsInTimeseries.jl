using TransitionsInTimeseries, Test

@testset "sliding ridge regression" begin
    n = 1001
    t = collect(1.0:n)
    x = copy(t)

    indicators = (mean, var)
    change_metric = RidgeRegressionSlope()

    config = WindowedIndicatorConfig(indicators, change_metric;
        width_ind = 100, stride_ind = 1,
        width_cha = 100, stride_cha = 1, whichtime = last,
    )

    res = estimate_transitions(config, x, t)

    # The trend of mean(windowview) is the stride for x=t
    meantrend_ground_truth = fill(1, length(res.t_change))
    # The trend of var(windowview) is 0 for x any affine function of t.
    vartrend_ground_truth = fill(0.0, length(res.t_change))
    @test isapprox(res.x_change[:, 1], meantrend_ground_truth, atol = 1e-9)
    @test isapprox(res.x_change[:, 2], vartrend_ground_truth, atol = 1e-9)

    # Virtually all results should have 0 significance versus the surrogates
    signif = SurrogatesSignificance(n = 100, tail = :both)
    pvalues = significant_transitions(res, signif)

    @test all(>(0.1), pvalues)

end