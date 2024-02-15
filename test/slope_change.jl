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

    tcross = res.t_change
    t1[end]

    @test tcross ≈ t1[end] atol = 1e-1
    @test b ≈ 1 atol = 1e-1
    @test d ≈ -1 atol = 1e-1
end

@testset "indicator" begin
    y = cumsum(x)
    config = SlopeChangeConfig(indicators = y -> y[2] - y[1], width_ind = 2)
    res = estimate_changes(config, y, t)
    @test tcross ≈ t1[end] atol = 1e-1
end