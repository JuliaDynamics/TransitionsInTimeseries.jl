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
    sol = solve(prob, EM(), dt=t[2]-t[1])

    return hcat(sol.(t)...)
end

@testset "ar1_whitenoise" begin
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
    # println([λ_est_vec, λ_est_mat, λ_est_gpu, λ])
    @test isapprox(λ_est_vec, λ_est_mat)
    @test isapprox(λ_est_gpu, λ_est_mat)
    @test isapprox(λ_est_gpu, λ, atol = T(1e-2))
end

function generate_affine_data(t::Vector{T}, nr::Int) where {T}
    Y = fill(T(NaN), (nr, length(t))) # init test data.
    W = 10 .* randn(2, nr)          # random parameters of affine function.
    for i in 1:nr
        Y[i, :] = W[1, i] .* t .+ W[2, i]
    end
    return Y, W
end

@testset "ridge_regression" begin
    T = Float32
    dt = T(1e-2) 
    t = collect(T(0):dt:T(10))
    nr = 10                         # number of tested regressions.
    Y, W = generate_affine_data(t, nr)
    Wcpu = ridge_regression(Y, t)
    Wgpu = Array( ridge_regression(CuArray(Y), t) )
    # println([W, Wcpu, Wgpu])
    @test sum(isapprox.(W, Wcpu, atol = T(1e-5))) == 2*nr
    @test sum(isapprox.(W, Wgpu, atol = T(1e-5))) == 2*nr
end