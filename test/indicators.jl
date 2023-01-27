using TransitionIndicators, Test
using TimeseriesSurrogates, Random

#################################################
# AR1 indicators
#################################################

# Check if AR1 regression parameter from a known AR1 process with white noise
# is successfully estimated.
@testset "ar1_whitenoise" begin
    θ = rand()
    x = AR1(;n_steps = 10000, x₀ = rand(), k = θ, rng = Random.default_rng())
    θ_est = ar1_whitenoise(x)
    @test isapprox(θ_est, θ, atol = 1e-3)
end
