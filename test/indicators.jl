using TransitionIndicators, Test, DifferentialEquations

#################################################
# AR1 indicators
#################################################

function f_ar1!(dx, x, p, t)
    dx[1] = p[1] * x[1]
end

function g_whitenoise!(dx, x, p, t)
    dx[1] = p[2]
end

function generate_ar1(x0::Vector{T}, p::Vector{T}, t::Vector{T}) where {T<:Real}
    tspan = extrema(t)
    prob = SDEProblem(f_ar1!, g_whitenoise!, x0, tspan, p)
    sol = solve(prob, EM(), dt=t[2]-t[1])
    return hcat(sol.(t)...)
end

# Check if AR1 regression parameter from a known AR1 process with white noise
# is successfully estimated.
@testset "ar1_whitenoise" begin
    T = Float32
    dt = T(1e-2)
    x0 = [T(1)]

    t = collect(T(0):dt:T(10))
    p = [-rand(T), T(0.1)]      # restoring rate λ (∈ [-1, 0] for stability), noise amplitude σ
    θ = exp(dt * p[1])          # θ is time discrete equivalent of λ

    x = vec(generate_ar1(x0, p, t))
    θ_est = ar1_whitenoise(x)

    @test isapprox(θ_est, θ, atol = T(1e-3))
end
