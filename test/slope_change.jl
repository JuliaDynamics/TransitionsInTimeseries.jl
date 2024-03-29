using TransitionsInTimeseries, Test, Statistics, Random

rng = Xoshiro(1234)
t1 = cumsum(rand(rng, 20))
t2 = cumsum(rand(rng, 20)) .+ t1[end]
b, d = 1.0, -1.0
x1 = b .* t1
x2 = x1[end] .+ d .* (t2 .- t1[end])
t = vcat(t1, t2)
x = vcat(x1, x2) .+ randn(rng, 40)./10

@testset "no indicator" begin
    config = SlopeChangeConfig()
    res = estimate_changes(config, x, t)
    a, b, c, d = res.fitparams

    tcross = first(res.t_change)

    @test tcross ≈ t1[end] atol = 1e-1
    @test b ≈ 1 atol = 1e-1
    @test d ≈ -1 atol = 1e-1
end

@testset "indicator" begin
    y = cumsum(x)
    config = SlopeChangeConfig(indicators = y -> y[2] - y[1], width_ind = 2, whichtime = last)
    res = estimate_changes(config, y, t)
    a, b, c, d = res.fitparams

    tcross = first(res.t_change)
    @test tcross ≈ t1[end] atol = 2e-1
    @test b ≈ 1 atol = 1e-1
    @test d ≈ -1 atol = 1e-1
end

@testset "significance" begin
    function fakedata(σ = 0.5)
        x1 = σ*randn(rng, 100)
        x2 = σ*randn(rng, 100) .+ range(0, 5; length = 100)
        # 1st: slope = 0, 2nd: slope = 0.05
        return vcat(x1, x2)
    end

    x = fakedata()
    for (flag, σ) in zip((true, false), (0.5, 5.0))
        x = fakedata(σ)
        res = estimate_changes(SlopeChangeConfig(), x)
        signif = SlopeChangeSignificance(; moe_slope = 0.1, moe_offset = 1.0,
            slope_diff = 0.03)  # according to only 0.05 slope difference in data
        out = significant_transitions(res, signif)
        @test first(out) == flag
    end
end
