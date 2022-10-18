using LinearAlgebra, CUDA
include("utils.jl")

#####################################################
#################### Windowing ######################
#####################################################
struct WindowingParams
    dt::Real
    Twndw::Real
    Tstrd::Real
    Nwndw::Int
    Nstrd::Int
end

"""
    get_windowing_params(Tvec::Vector{T})

Creates a WindowingParams struct out of a vector `[dt, Twndw, Tstrd]` with:
    - `dt` the time step of the time series to analyse.
    - `Twndw` the half-width of the window.
    - `Tstrd` the stride with which the windowing is applied.
"""
function get_windowing_params(Tvec::Vector{T}) where {T<:Real}
    N = get_step.(Tvec[2:end], Tvec[1])
    return WindowingParams(Tvec..., N...)
end

"""
    centered_wndw(X::AbstractArray, idx::Int, hw::Int)

Gets the center-windowed array of half-width `hw` at index `idx`
"""
centered_wndw(x::Vector{T}, idx::Int, hw::Int) where {T<:Real} = x[ (idx - hw):( idx + hw ) ]
centered_wndw(X::Matrix{T}, idx::Int, hw::Int) where {T<:Real} = X[ :, (idx - hw):( idx + hw ) ]
centered_wndw(X::CuArray{T, 2}, idx::Int, hw::Int) where {T<:Real} = X[ :, (idx - hw):( idx + hw ) ]

"""
    left_wndw(X::AbstractArray, idx::Int, hw::Int)

Gets the left-windowed array of half-width `hw` at index `idx`
"""
left_wndw(x::Vector{T}, idx::Int, hw::Int) where {T<:Real} = x[ (idx - 2*hw):idx ]
left_wndw(X::Matrix{T}, idx::Int, hw::Int) where {T<:Real} = X[ :, (idx - 2*hw):idx ]
left_wndw(X::CuArray{T, 2}, idx::Int, hw::Int) where {T<:Real} = X[ :, (idx - 2*hw):idx ]

"""
    right_wndw(X::AbstractArray, idx::Int, hw::Int)

Gets the right-windowed array of half-width `hw` at index `idx`
"""
right_wndw(X::Matrix{T}, idx::Int, hw::Int) where {T<:Real} = X[ :, idx:(idx + 2*hw) ]
right_wndw(x::Vector{T}, idx::Int, hw::Int) where {T<:Real} = x[ idx:(idx + 2*hw) ]
right_wndw(X::CuArray{T, 2}, idx::Int, hw::Int) where {T<:Real} = X[ :, idx:(idx + 2*hw) ]

# Multiple dispatch of the above-defined functions for internal use.
centered_wndw(n_wndw::Int, n_strd::Int, nt::Int) = (n_wndw+1):n_strd:(nt-n_wndw)
left_wndw(n_wndw::Int, n_strd::Int, nt::Int) = (2*n_wndw+1):n_strd:nt
right_wndw(n_wndw::Int, n_strd::Int, nt::Int) = 1:n_strd:(nt-2*n_wndw)


"""
    trim_wndw(X::AbstractArray, p::WindowingParams, wndw::Function)

Trims an array (of dim < 3) along its last dimension.
This is done in accordance with the windowing parameters `p` and the windowing function `wndw`.
"""
function trim_wndw(x::Vector{T}, p::WindowingParams, wndw::Function) where {T<:Real} 
    return x[ wndw(p.Nwndw, p.Nstrd, length(x)) ]
end

function trim_wndw(X::Matrix{T}, p::WindowingParams, wndw::Function) where {T<:Real}
    return X[ :, wndw(p.Nwndw, p.Nstrd, size(X, 2)) ]
end

#####################################################
################ Sliding estimators #################
#####################################################

"""
    slide_estimator(X::AbstractArray, p::WindowingParams, estimator::Function, wndw::Function)

Slides the computation of `estimator` over the last dimension of an array (dim < 3).
The type of windowing is specified by `wndw` and `p` (size and stride).
If X is a CuArray, the computation takes place on the GPU.
"""

function slide_estimator(
    x::Vector{T},
    p::WindowingParams,
    estimator::Function,
    wndw::Function,
) where {T<:Real}

    nt = length(x)
    strided_idx = wndw(p.Nwndw, p.Nstrd, nt)
    nidx = length(strided_idx)

    transition_indicator = fill(T(NaN), nidx)
    @inbounds for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        transition_indicator[j1] = estimator( wndw( x, j2, p.Nwndw ) )
    end
    return transition_indicator
end

function slide_estimator(
    X::Matrix{T},
    p::WindowingParams,
    estimator::Function,
    wndw::Function,
) where {T<:Real}

    nl, nt = size(X)
    strided_idx = wndw(p.Nwndw, p.Nstrd, nt)
    nidx = length(strided_idx)

    transition_indicator = fill(T(NaN), nl, nidx)
    @inbounds for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        transition_indicator[:, j1] = estimator( wndw( X, j2, p.Nwndw ) )
    end
    return transition_indicator
end

function slide_estimator(
    X::CuArray{T, 2},
    p::WindowingParams,
    estimator::Function,
    wndw::Function,
) where {T<:Real}

    nl, nt = size(X)
    strided_idx = wndw(p.Nwndw, p.Nstrd, nt)
    nidx = length(strided_idx)

    transition_indicator = fill(T(NaN), nl, nidx)
    @inbounds for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        transition_indicator[:, j1] = Array( estimator( wndw( X, j2, p.Nwndw ) ) )
    end
    return CuArray( transition_indicator )
end

#####################################################
#################### Smoothing ######################
#####################################################

"""
    gettrend_rollmean(x::Vector{Real}, p::WindowingParams, wndw::Function)

Computes the rolling mean of the vector x.
This is done according to the windowing parameters `p` and the function `wndw`.
"""
function gettrend_rollmean(x::Vector{T}, p::WindowingParams, wndw::Function) where {T<:Real}
    return slide_estimator(x, p, mean, wndw)
end

"""
    gettrend_rollkernel(x::Union{Vector{Real}, Matrix{Real}}, p::WindowingParams, wndw::Function, kernel::Function)

Computes the kernel convolution of x over its last dimension.
This is done according to the windowing parameters `p` and the function `wndw`.
Functions available for `kernel`:
    - uniform_kernel
    - triangular_kernel
    - parabolic_kernel
    - biweight_kernel
    - triweight_kernel
    - tricube_kernel
    - gaussian_kernel
    - cosine_kernel
    - logistic_kernel
    - sigmoid_kernel
"""
function gettrend_rollkernel(
    x::Vector{T},
    p::WindowingParams,
    wndw::Function,
    kernel::Function,
) where {T<:Real}
    
    kernel = scaled_kernel(T, p, kernel)
    nt = length(x)
    strided_idx = wndw(p.Nwndw, p.Nstrd, nt)
    nidx = length(strided_idx)

    filtered_signal = fill(T(NaN), nidx)
    @inbounds for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        filtered_signal[j1] = wndw( x, j2, p.Nwndw )' .* kernel
    end
    return filtered_signal
end

function gettrend_rollkernel(
    X::Matrix{T},
    p::WindowingParams,
    wndw::Function,
    kernel::Function,
) where {T<:Real}

    kernel = scaled_kernel(T, p, kernel)
    nl, nt = size(X)
    strided_idx = wndw(p.Nwndw, p.Nstrd, nt)
    nidx = length(strided_idx)

    filtered_signal = fill(T(NaN), nl, nidx)
    @inbounds for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        filtered_signal[:, j1] = wndw( X, j2, p.Nwndw )' .* kernel
    end
    return filtered_signal
end

"""
    detrend(x::Vector{T}, x_trend::Vector{T})
"""
function detrend(x::Vector{T}, x_trend::Vector{T}) where {T<:Real}
    return x - x_trend
end

# TODO implement further ways to detrend the signal such as the one below.
function gettrend_loess(x::Vector{T}, σ::Int) where {T<:Real}
end

function gettrend_dfa(x::Vector{T}, σ::Int) where {T<:Real}
end

function gettrend_emd(x::Vector{T}, σ::Int) where {T<:Real}
end

