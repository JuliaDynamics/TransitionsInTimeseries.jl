using Downloads, DelimitedFiles
using TransitionIndicators, Test, Random, TimeseriesSurrogates

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
    indconfig = IndicatorsConfig(t, [LowfreqPowerSpectrum], width = nt)

    @test indconfig.indicators[1](y_lofreq[1:indconfig.width]) > 0.95
    @test indconfig.indicators[1](y_hifreq[1:indconfig.width]) < 0.05
end

# TODO: improve this test. For now, simply based on what is observed in:
# https://github.com/JuliaDynamics/NonlinearDynamicsComplexSystemsCourse/blob/main/notebooks/nonlinear_causal_timeseries_analysis.ipynb
@testset "perment" begin
    file = Downloads.download("https://raw.githubusercontent.com/JuliaDynamics/" * 
        "NonlinearDynamicsTextbook/master/exercise_data/7.csv")
    x = vec(DelimitedFiles.readdlm(file))
    t = collect(eachindex(x))

    indconfig = IndicatorsConfig(t, [PermutationEntropy], width = 50)
    sigconfig = SignificanceConfig(indconfig, [RidgeRegressionSlope], width = 50, stride = 1, n_surrogates = 2)
    res = indicators_analysis(t, x, indconfig, sigconfig)
    @test mean(res.x_indicator[7000:7500]) < mean(res.x_indicator[5500:6000])
end