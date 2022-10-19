using Test, TransitionIdentifiers
using DifferentialEquations, CUDA, CairoMakie

function f_linear!(dx, x, p, t)
    dx[1] = p[1] * x[1]
end

function g_whitenoise!(dx, x, p, t)
    dx[1] = p[2]
end

function generate_ar1(x0::Vector{T}, p::Vector{T}, t::Vector{T}) where {T}
    tspan = extrema(t)
    prob = SDEProblem(f_linear!, g_whitenoise!, x0, tspan, p)
    sol = solve(prob)
    return hcat(sol.(t)...)
end

@testset begin
    T = Float32
    dt = T(1e-2) 
    t = collect(T(0):dt:T(10))
    p = [-rand(T), T(0.1)]      # λ, σ
    x0 = [T(1)]

    x_mat = generate_ar1(x0, p, t)
    x_vec = vec(x_mat)
    x_gpu = CuArray(x_mat)
    λ_est_vec = ar1_whitenoise(x_vec)
    λ_est_mat = ar1_whitenoise(x_mat)[1]
    λ_est_gpu = Array( ar1_whitenoise(x_gpu) )[1]
    λ = exp(dt * p[1])
    λ_list = [λ_est_vec, λ_est_mat, λ_est_gpu, λ]
    println(λ_list)
    @test isapprox(λ_est_vec, λ_est_mat)
    @test isapprox(λ_est_gpu, λ_est_mat)
    @test isapprox(λ_est_gpu, λ, atol = T(1e-2))
end

@testset begin
    dt = 1f-2
    t = collect(0f0:dt:10f0)
    tspan = extrema(t)
    a = rand()
end