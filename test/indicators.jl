using TransitionsInTimeseries, Test, Random, TimeseriesSurrogates, Distributions

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

# Test kolmogorov_smirnov by sampling different distributions
@testset "kolmogorov_smirnov" begin
    n = 1000

    distributions = [Uniform(), Normal(), Binomial()]
    for (i, d1) in enumerate(distributions)
        for (j, d2) in enumerate(distributions)
            x = vcat(rand(d1, n), rand(d2, n))
            if i == j
                @test kolmogorov_smirnov(x) > 0.1
            else
                @test kolmogorov_smirnov(x) < 1e-8
            end
        end
    end
end