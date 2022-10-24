#####################################################
# Fourier surrogates
#####################################################

# Shift each frequency content by a random phase.
function generate_fourier_surrogate(x::Vector{T}) where {T<:Real}
    F = rfft(x)
    return irfft(F .* exp.(2 * π * im .* rand(length(F))), length(x))
end

"""
    generate_stacked_fourier_surrogates(X::AbstractArray, ns::Int)

Generate ns Fourier surrogates for each line of X.
If X is a CuArray, computation takes place in GPU.
Aformentionned speedup particularly significant if X big or ns large.
"""
function generate_stacked_fourier_surrogates(x::Vector{T}, ns::Int) where {T<:Real}
    S = zeros(T, ns, length(x))
    for i in axes(S, 1)
        S[i, :] = generate_fourier_surrogate(x)
    end
    return S
end

function generate_stacked_fourier_surrogates(X::Matrix{T}, ns::Int) where {T<:Real}
    nx, nt = size(X)
    F = repeat(rfft(X, 2), inner = (ns, 1))
    stacked_surrogates = irfft(
        F .* exp.(T(2 * π) * im .* rand(T, nx * ns, size(F, 2))),
        nt,
        2,
    )
    return stacked_surrogates
end

function generate_stacked_fourier_surrogates(X::CuArray{T,2}, ns::Int) where {T<:Real}
    nx, nt = size(X)
    F = repeat(CUDA.CUFFT.rfft(X, 2), inner = (ns, 1))
    stacked_surrogates = CUDA.CUFFT.irfft(
        F .* exp.(T(2 * π) * im .* CUDA.rand(T, nx * ns, size(F, 2))),
        nt,
        2,
    )
    return stacked_surrogates
end

#####################################################
# Ridge regression
#####################################################

"""
    ridge_regression(y::AbstractArray, t::Vector; lambda=0::Real)

Ridge regression with regularization parameter lambda.
Set lambda=0 to recover linear regression.
"""
function ridge_regression(y::Vector{T}, t::Vector{T}; lambda = 0::Real) where {T<:Real}
    t = t .- t[1]
    T_bias_ext = hcat(t, ones(length(t)))'
    return inv(T_bias_ext * T_bias_ext' + lambda .* I(2)) * T_bias_ext * y
end

function ridge_regression(Y::Matrix{T}, t::Vector{T}; lambda = 0::Real) where {T<:Real}
    t = t .- t[1]
    T_bias_ext = hcat(t, ones(length(t)))'
    return inv(T_bias_ext * T_bias_ext' + lambda .* I(2)) * T_bias_ext * transpose(Y)
end

# GPU version of ridge regression. Y must be nt x ns.
function ridge_regression(Y::CuArray{T,2}, t::Vector{T}; lambda = 0::Real) where {T<:Real}
    t = t .- t[1]
    T_bias_ext = hcat(t, ones(length(t)))'
    return CuArray(inv(T_bias_ext * T_bias_ext' + lambda .* I(2))) *
           CuArray(T_bias_ext) *
           permutedims(Y)
end

"""
    ridge_regression_slope(y::AbstractArray, t::Vector; lambda=0::Real)

Get the slope of underlying ridge regression.
"""
function ridge_regression_slope(
    Y::Matrix{T},
    t::Vector{T};
    lambda = 0::Real,
) where {T<:Real}
    return ridge_regression(Y, t; lambda)[1, :]
end

function ridge_regression_slope(
    Y::CuArray{T,2},
    t::Vector{T};
    lambda = 0::Real,
) where {T<:Real}
    return ridge_regression(Y, t; lambda)[1, :]
end

#####################################################
# Kendall tau
#####################################################

function kendall_tau(y::Vector{T}, t::Vector{T}) where {T<:Real}
    return StatsBase.corkendall(t, y)
end

# This implementation works fine for even time spacing
# because the time vector does not need to be shifted.
# TODO: make more general!
function kendall_tau(Y::Matrix{T}, t::Vector{T}) where {T<:Real}
    return mapslices(x -> StatsBase.corkendall(t, x), Y, dims = 2)
end

# Scaling with mean (or other function) might be a better metric than
# just the trend! Typically: AR1 needs to be close to 1.
function scaled_kendall_tau(Y::Matrix{T}, t::Vector{T}; scaling::Function=StatsBase.mean) where {T<:Real}
    return mapslices(x -> StatsBase.corkendall(t, x) * scaling(x), Y, dims = 2)
end

#####################################################
# Percentile significance 
#####################################################
"""
    percentile_significance(ref_stat::AbstractArray, sur_stat::Matrix{T}, ns::Int, nx::Int)

Get percentile significance of reference (original time series) versus surrogate statistics.
"""
function percentile_significance(
    ref_stat::Matrix{T},
    sur_stat::Matrix{T},
    ns::Int,
    nx::Int,
) where {T<:Real}

    rep_ref = repeat(ref_stat, inner = (ns, 1))
    return (kron(I(nx), ones(T, ns)') * (sur_stat .- rep_ref .< T(0))) ./ T(ns)
end

function percentile_significance(
    ref_stat::CuArray{T,2},
    sur_stat::CuArray{T,2},
    ns::Int,
    nx::Int,
) where {T<:Real}

    rep_ref = repeat(ref_stat, inner = (ns, 1))
    return (CuArray(kron(I(nx), ones(T, ns)')) * (sur_stat .- rep_ref .< T(0))) ./ T(ns)
end

#####################################################
# Slide trend estimations
#####################################################
"""
    slide_idtrend(
        X::CuArray{T, 2},
        t::Vector{T},
        p::WindowingParams,
        estimator::Function,
        wndw::Function;
        kwargs...,
    )

Slide trend estimation over computed identifier time series `X`. 
"""
function slide_idtrend(
    X::CuArray{T,2},
    t::Vector{T},
    p::WindowingParams,
    estimator::Function,
    wndw::Function;
    kwargs...,
) where {T<:Real}

    nl, nt = size(X)
    strided_idx = wndw(p.Nwndw, p.Nstrd, nt)
    nidx = length(strided_idx)

    idtrend = fill(T(NaN), nl, nidx)
    for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        idtrend[:, j1] =
            Array(estimator(wndw(X, j2, p.Nwndw), wndw(t, j2, p.Nwndw); kwargs...))
    end
    return CuArray(idtrend)
end

function slide_idtrend(
    X::Matrix{T},
    t::Vector{T},
    p::WindowingParams,
    estimator::Function,
    wndw::Function;
    kwargs...,
) where {T<:Real}

    nl, nt = size(X)
    strided_idx = wndw(p.Nwndw, p.Nstrd, nt)
    nidx = length(strided_idx)

    idtrend = fill(T(NaN), nl, nidx)
    for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        idtrend[:, j1] = estimator(wndw(X, j2, p.Nwndw), wndw(t, j2, p.Nwndw); kwargs...)
    end
    return idtrend
end

#####################################################
# Predict transition
#####################################################

# TODO insert tolerance wrt lag
function predict_transition(
    indicator_list::Vector{Matrix{T}};
    plevel=T(0.95),
    nindicators::Int=length(indicator_list)
) where {T}
    P = stack_indicators(indicator_list)
    S = sum_significant_indicators(P, plevel)
    return threshold_indicator_significance(S, nindicators)
end

function stack_indicators(indicator_list::Vector{Matrix{T}}) where {T}
    return cat(indicator_list..., dims=3)
end

function sum_percentiles(P::Array{T, 3}) where {T}
    return reduce( +, P, dims=3 )
end

function sum_significant_indicators(P::Array{T, 3}, plevel::T) where {T}
    return reduce( +, P .> plevel, dims=3 )
end

function threshold_indicator_significance(S::Array{T, 3}, nindicators::Int) where {T}
    return (S .>= nindicators)[:,:,1]
end