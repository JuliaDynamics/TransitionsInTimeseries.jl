using TransitionsInTimeseries, Test, Statistics, Random, Distributions

@testset "missmatch in length" begin
    indicators = (mean, var)
    change_metric = RidgeRegressionSlope()
    @test_throws ArgumentError SlidingWindowConfig(indicators, change_metric)
    @test_throws ArgumentError SegmentedWindowConfig(indicators, change_metric, [1], [100])
end

@testset "sliding ridge regression" begin
    n = 10001
    t = collect(1.0:n)
    x = copy(t)
    w = 100
    s = 1
    sqsum(x) = sum(x.^2)    # used only for QuantileSignificance

    indicators = (mean, var, sqsum)
    change_metric = (RidgeRegressionSlope(), RidgeRegressionSlope(), RidgeRegressionSlope())

    config = SlidingWindowConfig(indicators, change_metric;
        width_ind = w, stride_ind = s,
        width_cha = w, stride_cha = s, whichtime = last,
    )

    res = estimate_changes(config, x, t)
    ni, nc = length(res.t_indicator), length(res.t_change)

    # Mean of identity is identity with offset; its trend is the stride for x=t
    mean_ground_truth = range((w+1)/2, step = s, length = length(t) - w + 1)
    meantrend_ground_truth = fill(1, nc)

    # Var of identity is constant and trend is 0.
    var_groundtruth = fill(var(1:w), ni)
    vartrend_ground_truth = fill(0.0, nc)

    @test isapprox(res.x_indicator[:, 1], mean_ground_truth, atol = 1e-9)
    @test isapprox(res.x_change[:, 1], meantrend_ground_truth, atol = 1e-9)
    @test isapprox(res.x_indicator[:, 2], var_groundtruth, atol = 1e-9)
    @test isapprox(res.x_change[:, 2], vartrend_ground_truth, atol = 1e-9)

    # Virtually all results should have 0 significance versus the surrogates
    # Compared to segmnented window: increase of mean not significant over such small windows
    signif = SurrogatesSignificance(n = 100, tail = [:both, :both, :both], p = 0.1)
    flags = significant_transitions(res, signif)

    @test all(.!flags)

    # Trend in mean = 1, trend in var = 0. Check separation.
    signif = ThresholdSignificance(0.5)
    flags = significant_transitions(res, signif)
    @test flags[:, 1] == fill(true, nc)
    @test flags[:, 2] == fill(false, nc)

    # Since mean and var have constant trend, they are not suited
    # to test the quantile significance, which is based on the spread of x_change.
    # We therefore use the result of sqsum, which should display 10% of positives
    # with a an absolute tolerance of 2 for the default parameters of QuantileSignificance
    signif = QuantileSignificance()
    flags = significant_transitions(res, signif)
    min_pos_amount = floor(0.1 * nc)
    @test min_pos_amount <= sum(flags[:, 3]) <= min_pos_amount + 2

    # We should roughly have 5% positives for a normal distributions and 2σ signif.
    view(res.x_change, :, 3) .= rand(Normal(), size(res.x_change, 1))
    signif = SigmaSignificance(factor = 2)
    flags = significant_transitions(res, signif)
    min_pos_amount, max_pos_amount = floor(0.04 * nc), floor(0.06 * nc)
    @test min_pos_amount <= sum(flags[:, 3]) <= max_pos_amount

end

@testset "segmented ridge regression" begin
    n = 1001
    t = collect(1.0:n)
    x = copy(t)
    w = 100
    s = 1

    indicators = (mean, var)
    change_metric = (RidgeRegressionSlope(), RidgeRegressionSlope())
    config = SegmentedWindowConfig(indicators, change_metric, [t[1]], [t[end]];
        width_ind = w, stride_ind = s, min_width_cha = 30, whichtime = last,
    )

    res = estimate_changes(config, x, t)

    # Mean of identity is identity with offset; its trend is the stride for x=t
    mean_ground_truth = range((w+1)/2, step = s, length = length(t) - w + 1)
    meantrend_ground_truth = 1

    # Var of identity is constant and trend is 0.
    var_groundtruth = fill(var(1:w), length(res.t_indicator[1]))
    vartrend_ground_truth = 0

    @test isapprox(res.x_indicator[1][:, 1], mean_ground_truth, atol = 1e-9)
    @test isapprox(res.x_change[1, 1], meantrend_ground_truth, atol = 1e-9)
    @test isapprox(res.x_indicator[1][:, 2], var_groundtruth, atol = 1e-9)
    @test isapprox(res.x_change[1, 2], vartrend_ground_truth, atol = 1e-9)

    # Virtually all results should have 0 significance versus the surrogates
    signif = SurrogatesSignificance(n = 1000, tail = [:both, :both], p = 0.1)
    flags = significant_transitions(res, signif)

    # Increase in mean is significant, increase in var is not
    @test flags[1, 1] == true
    @test flags[1, 2] == false
end

@testset "nothing indicator" begin
    rng = Xoshiro(123)

    N = 1000 # the statistic is independent of `N` for large enough `N`!
    x = randn(rng, N)
    y = 1.8randn(rng, N) .+ 4.0
    z = vcat(x, y)

    config = SlidingWindowConfig(nothing, (difference_of_means, difference_of_maxes);
        width_cha = 100, stride_cha = 50)

    res = estimate_changes(config, z)

    c = res.x_change[:, 1]

    @test maximum(c) > 1
    @test count(>(1), c) == 1

    # surrogates
    signif = SurrogatesSignificance(surromethod = RandomShuffle(), n = 1000, tail = [:right, :right], p = 0.05)
    flags = significant_transitions(res, signif)

    @test count(flags) == 2

end