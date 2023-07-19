using TransitionsInTimeseries, Test, Random, TimeseriesSurrogates

# Check if AR1 regression parameter from a known AR1 process with white noise
# is successfully estimated.
@testset "ar1_whitenoise" begin
    θ = rand()
    x = AR1(1_000_000, rand(), θ, Random.default_rng())
    θ_est = ar1_whitenoise(x)
    @test isapprox(θ_est, θ, atol = 5e-3)
end

@testset "lfps" begin
    t = 0:0.1:100
    nt = length(t)
    y_lofreq = fill(10 * rand(), nt)
    y_hifreq = sin.(t .* 100)
    metric = LowfreqPowerSpectrum()
    metricc = precompute(metric, ones(nt))

    @test metricc(y_lofreq[1:nt]) > 0.95
    @test metricc(y_hifreq[1:nt]) < 0.05
end

 @testset "perment" begin
    x = 1:1000
    pe = permutation_entropy()
    # Sliding window over this x gives the same entropy values, all of which are zero
    res = windowmap(pe, x; width = 10)
    @test res == zeros(length(res))
end