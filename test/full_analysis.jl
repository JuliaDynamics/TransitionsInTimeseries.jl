using TransitionsInTimeseries, Test, Statistics, Random

@testset "sliding ridge regression" begin
    n = 1001
    t = collect(1.0:n)
    x = copy(t)

    indicators = (mean, var)
    change_metric = (RidgeRegressionSlope(), RidgeRegressionSlope())

    config = SlidingWindowConfig(indicators, change_metric;
        width_ind = 100, stride_ind = 1,
        width_cha = 100, stride_cha = 1, whichtime = last,
    )

    res = estimate_indicator_changes(config, x, t)

    # The trend of mean(windowview) is the stride for x=t
    meantrend_ground_truth = fill(1, length(res.t_change))
    # The trend of var(windowview) is 0 for x any affine function of t.
    vartrend_ground_truth = fill(0.0, length(res.t_change))
    @test isapprox(res.x_change[:, 1], meantrend_ground_truth, atol = 1e-9)
    @test isapprox(res.x_change[:, 2], vartrend_ground_truth, atol = 1e-9)

    # Virtually all results should have 0 significance versus the surrogates
    signif = SurrogatesSignificance(n = 100, tail = :both, p = 0.1)
    flags = significant_transitions(res, signif)

    @test all(.!flags)

    @testset "missmatch in length" begin
        indicators = (mean, var)
        change_metric = RidgeRegressionSlope()
        @test_throws ArgumentError SlidingWindowConfig(indicators, change_metric)
    end

end

@testset "nothing indicator" begin
    rng = Xoshiro(123)

    N = 1000 # the statistic is independent of `N` for large enough `N`!
    x = randn(rng, N)
    y = 1.8randn(rng, N) .+ 4.0
    z = vcat(x, y)

    config = SlidingWindowConfig(nothing, (difference_of_means, difference_of_maxes);
        width_cha = 100, stride_cha = 50)

    res = estimate_indicator_changes(config, z)

    c = res.x_change[:, 1]

    @test maximum(c) > 1
    @test count(>(1), c) == 1

    # surrogates
    signif = SurrogatesSignificance(surromethod = RandomShuffle(), n = 1000, tail = :right, p = 0.05)
    flags = significant_transitions(res, signif)

    @test count(flags) == 2

end