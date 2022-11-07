#=
TODO: test statistical moments (wrt. StatsBase)
=#

using Revise
push!(LOAD_PATH,"/home/jan/.julia/dev/TransitionIdentifiers")
using TransitionIdentifiers
using Test
using DifferentialEquations
using CUDA
using StatsBase

###############################
# Correct computation of statistical moments
###############################

@testset "skewness" begin
    T = Float32
    n = 10000
    x = rand(T, n)
    X = reshape(x, (1,n))
    Xgpu = CuArray(X)

    skw0 = StatsBase.skewness(x)
    skw1 = skw(X)[1, 1]
    skw2 = skw(Xgpu)[1, 1]
    @test isapprox(skw0, skw1, atol = T(1e-5))
    @test isapprox(skw0, skw2, atol = T(1e-5))
end

@testset "kurtosis" begin
    T = Float32
    n = 10000
    x = rand(T, n)
    X = reshape(x, (1,n))
    Xgpu = CuArray(X)

    krt0 = StatsBase.kurtosis(x)
    krt1 = krt(X)[1, 1]
    krt2 = krt(Xgpu)[1, 1]
    @test isapprox(krt0, krt1, atol = T(1e-2))
    @test isapprox(krt0, krt2, atol = T(1e-2))
end

###############################
# Analytic AR1 regression
###############################

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
    p = [-rand(T), T(0.1)]      # λ (make it random ), σ
    x0 = [T(1)]
    θ = exp(dt * p[1])          # θ is time discrete equivalent of λ

    x_mat = generate_ar1(x0, p, t)
    x_vec = vec(x_mat)
    x_gpu = CuArray(x_mat)
    θ_est_vec = ar1_whitenoise(x_vec)
    θ_est_mat = ar1_whitenoise(x_mat)[1]
    θ_est_gpu = Array( ar1_whitenoise(x_gpu) )[1]
    # θ_est_msk = Array( masked_ar1_whitenoise(xgpu) )
    # TODO: add the accelerated AR1

    @test isapprox(θ_est_vec, θ, atol = T(1e-3))
    @test isapprox(θ_est_mat, θ, atol = T(1e-3))
    @test isapprox(θ_est_gpu, θ, atol = T(1e-3))
    # @test isapprox(θ_est_msk, θ, atol = T(1e-3))
end

###############################
# Linear regression
###############################

function generate_affine_data(t::Vector{T}, nr::Int) where {T}
    Y = fill(T(NaN), (nr, length(t)))   # init test data.
    W = 10 .* randn(2, nr)              # random parameters of affine function.
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
    @test sum(isapprox.(W, Wcpu, atol = T(1e-5))) == length(W)
    @test sum(isapprox.(W, Wgpu, atol = T(1e-5))) == length(W)
end

###############################
# Percentile computation
###############################

@testset "significance" begin
    T = Float32
    nx, ns, nt = 2, 100, 10         # n° of variables, surrogates per variable, time steps per variable.
    sur_stat = repeat( collect(1:ns) * ones(T, nt)', outer=(2,1) )
    ref_stat = round.( rand(T, nx, nt) * 100 )
    significance = percentile_significance(ref_stat, sur_stat, ns, nx)
    @test significance .* 100 .+ 1 == ref_stat
end