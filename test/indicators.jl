using TransitionIndicators, Test, Random, TimeseriesSurrogates

# Check if AR1 regression parameter from a known AR1 process with white noise
# is successfully estimated.
@testset "ar1_whitenoise" begin
    θ = rand()
    x = AR1(100_000, rand(), θ) #, Random.default_rng())
    θ_est = ar1_whitenoise(x)
    @test isapprox(θ_est, θ, atol = 5e-3)
end

@testset "lfps" begin
    t = 0:0.1:100
    nt = length(t)
    y_lofreq = fill(10 * rand(), nt)
    y_hifreq = sin.(t .* 100)
    iparams = IndicatorsParams(q_lofreq = 0.1)
    indconfig = IndicatorsConfig(t, [LowfreqPowerSpectrum], width = nt)

    @test indconfig.indicators[1](y_lofreq[1:indconfig.width]) > 0.95
    @test indconfig.indicators[1](y_hifreq[1:indconfig.width]) < 0.05
end