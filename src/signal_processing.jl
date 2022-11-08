#=
For the recognition of transition, following routines are necessary:
- Windowing
- Sliding functions over time series
- Smoothing and detrending
- Masking (only for high-performance on GPU)

They are implemented here.
=#

#####################################################
# Windowing
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
centered_wndw(x::Vector{T}, idx::Int, hw::Int) where {T<:Real} = x[(idx-hw):(idx+hw)]
centered_wndw(X::Matrix{T}, idx::Int, hw::Int) where {T<:Real} = X[:, (idx-hw):(idx+hw)]
centered_wndw(X::CuArray{T,2}, idx::Int, hw::Int) where {T<:Real} = X[:, (idx-hw):(idx+hw)]

"""

    left_wndw(X::AbstractArray, idx::Int, hw::Int)

Gets the left-windowed array of half-width `hw` at index `idx`
"""
left_wndw(x::Vector{T}, idx::Int, hw::Int) where {T<:Real} = x[(idx-2*hw):idx]
left_wndw(X::Matrix{T}, idx::Int, hw::Int) where {T<:Real} = X[:, (idx-2*hw):idx]
left_wndw(X::CuArray{T,2}, idx::Int, hw::Int) where {T<:Real} = X[:, (idx-2*hw):idx]

"""

    right_wndw(X::AbstractArray, idx::Int, hw::Int)

Gets the right-windowed array of half-width `hw` at index `idx`
"""
right_wndw(X::Matrix{T}, idx::Int, hw::Int) where {T<:Real} = X[:, idx:(idx+2*hw)]
right_wndw(x::Vector{T}, idx::Int, hw::Int) where {T<:Real} = x[idx:(idx+2*hw)]
right_wndw(X::CuArray{T,2}, idx::Int, hw::Int) where {T<:Real} = X[:, idx:(idx+2*hw)]

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
    return x[wndw(p.Nwndw, p.Nstrd, length(x))]
end

function trim_wndw(X::Matrix{T}, p::WindowingParams, wndw::Function) where {T<:Real}
    return X[:, wndw(p.Nwndw, p.Nstrd, size(X, 2))]
end

#####################################################
# Sliding estimators
#####################################################

"""

    slide(
        X::A,
        p::WindowingParams,
        estimator::Function,
        wndw::Function;
        kwargs...,
    ) where {A<:Union{Matrix{T}, CuArray{T, 2}}} where {T<:Real}

Function sliding `estimator` over time axis of `X` with time `t`, windowing parameters `p` and window type `wndw`.
"""
function slide(
    X::A,
    t::Vector{T},
    p::WindowingParams,
    estimator::Function,
    wndw::Function;
    kwargs...,
) where {A<:Union{Matrix{T}, CuArray{T, 2}}} where {T<:Real}

    nl, nt = size(X)
    strided_idx = wndw(p.Nwndw, p.Nstrd, nt)
    nidx = length(strided_idx)

    # Initialize result of sliding estimator over X.
    Y = fill(T(NaN), nl, nidx)
    @inbounds for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        Y[:, j1] =
            Array(estimator(wndw(X, j2, p.Nwndw), wndw(t, j2, p.Nwndw); kwargs...))
    end
    return A(Y)
end

function slide(
    x::A,
    t::Vector{T},
    p::WindowingParams,
    estimator::Function,
    wndw::Function;
    kwargs...,
) where {A<:Union{Vector{T}, CuArray{T, 1}}} where {T<:Real}

    X = reshape(x, (1, length(x)))
    return slide(X, t, p, estimator, wndw)
end

"""

    function grow_window(
        X::Matrix{T},
        t::Vector{T},
        p::WindowingParams,
        estimator::Function,
        wndw::Function;
        kwargs...,
    ) where {T<:Real}
"""
function grow_window(
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

    Y = fill(T(NaN), nl, nidx)
    for j1 in eachindex(strided_idx)
        j2 = strided_idx[j1]
        Y[:, j1] = estimator(X[:, 1:j2], t[1:j2]; kwargs...)
    end
    return Y
end
# TODO: test

#####################################################
# Smoothing
#####################################################

"""

    gettrend_rollmean(x::Vector{Real}, p::WindowingParams, wndw::Function)

Computes the rolling mean of the vector x.
This is done according to the windowing parameters `p` and the function `wndw`.
"""
function gettrend_rollmean(x::Vector{T}, p::WindowingParams, wndw::Function) where {T<:Real}
    return slide_estimator(x, p, StatsBase.mean, wndw)
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
        filtered_signal[j1] = wndw(x, j2, p.Nwndw)' * kernel
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
        filtered_signal[:, j1] = wndw(X, j2, p.Nwndw)' * kernel
    end
    return filtered_signal
end

function get_trend(
    X::Matrix{T},
    p::WindowingParams,
    wndw::Function,
    kernel::Function,
) where {T}
    return mapslices(x_ -> gettrend_rollkernel(x_, p, wndw, kernel), X; dims = 2)
end
# TODO replace mapslices by matrix vector multiplication (which can also be performed on GPU)

"""

    detrend(x::Vector{T}, x_trend::Vector{T})
"""
function detrend(x::Vector{T}, x_trend::Vector{T}) where {T<:Real}
    return x - x_trend
end

# TODO implement further ways to detrend the signal such as the one below.
function gettrend_loess(x::Vector{T}, σ::Int) where {T<:Real} end

function gettrend_dfa(x::Vector{T}, σ::Int) where {T<:Real} end

function gettrend_emd(x::Vector{T}, σ::Int) where {T<:Real} end

#####################################################
# Masking
#####################################################

"""

    window_mask(nt::Int, p::WindowingParams, wndw::Function)

Get matrix of type:
1 0 0
1 1 0
0 1 1
0 0 1

Windowing differentiation performed wrongly so far --> just append zero matrix for correct shift.
"""
function window_mask(T::Type, nt::Int, p::WindowingParams, type::String)
    if type == "common"
        idx = -p.Nwndw:p.Nwndw
    elseif type == "ar1"
        idx = -p.Nwndw+1:p.Nwndw
    end
    return spdiagm([i => ones(T, nt - abs(i)) for i in idx]...)
end

"""

    strided_window_mask(nt::Int, p::WindowingParams, wndw::Function)

Get matrix of type:
1 0
1 0
0 1
0 1

For AR1 computation, window needs to be adapted --> use e.g. ar1_left_wndw instead of left_wndw
"""
function strided_window_mask(
    T::Type,
    nt::Int,
    p::WindowingParams,
    type::String,
    wndw::Function,
)
    return window_mask(T, nt, p, type)[:, wndw(p.Nwndw, p.Nstrd, nt)]
end

# Multiple dispatch of windowing functions for mask generation.
function centered_wndw(M::SparseMatrixCSC{T,Int}, p::WindowingParams) where {T}
    n1, n2 = size(M)
    Z = spzeros(T, (n1, p.Nwndw))
    return hcat(Z, M[:, p.Nwndw+1:end-p.Nwndw], Z)
end

function left_wndw(M::SparseMatrixCSC{T,Int}, p::WindowingParams) where {T}
    n1, n2 = size(M)
    Z = spzeros(T, (n1, 2 * p.Nwndw))
    return hcat(Z, M[:, p.Nwndw+1:end-p.Nwndw])
end

function right_wndw(M::SparseMatrixCSC{T,Int}, p::WindowingParams) where {T}
    n1, n2 = size(M)
    Z = spzeros(T, (n1, 2 * p.Nwndw))
    return hcat(M[:, p.Nwndw+1:end-p.Nwndw], Z)
end
