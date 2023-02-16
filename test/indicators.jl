using TransitionIndicators, Test, Random, TimeseriesSurrogates

#################################################
# AR1 indicators
#################################################

# Check if AR1 regression parameter from a known AR1 process with white noise
# is successfully estimated.
@testset "ar1_whitenoise" begin
    θ = rand()
    x = AR1(1_000_000, rand(), θ, Random.default_rng())
    θ_est = ar1_whitenoise(x)
    @test isapprox(θ_est, θ, atol = 5e-3)
end
