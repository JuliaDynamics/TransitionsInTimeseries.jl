using TransitionIndicators, Test, Random

#################################################
# AR1 indicators
#################################################

# Check if AR1 regression parameter from a known AR1 process with white noise
# is successfully estimated.
@testset "ar1_whitenoise" begin
    θ = rand()
    x = AR1(;n_steps = 10_000, x₀ = rand(), k = θ)
    θ_est = ar1_whitenoise(x)
    @test isapprox(θ_est, θ, atol = 1e-2)
end
