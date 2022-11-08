#=
Transition indicators can be classified in two families:
- temporal
- spatial

The temporal indicators can be classified depending on the undlerying sampling:
- even
- uneven

... or depending on the type of computation they imply:
- statistical moments
- regression
- frequency analysis

For the latter classification, the most common metrics are implemented here.
TODO: implement temporal indicators on uneven sampling.
TODO: implement spatial indicators.

Finally, regression models rely on a noise assumption. The latter can be:
- white
- an AR1 process
- a higher-order process which is jointly estimated during the regression
TODO: implement beyond white-noise assumption.
=#

#####################################################
# Statistical moments
#####################################################

"""

    mean_lastdim(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the mean of an array over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""
function mean_lastdim(X::A) where {A<:Union{Array{T},CuArray{T}}} where {T<:Real}
    lastdim = length(size(X))
    return reduce(+, X, dims = lastdim) ./ T(size(X, lastdim))
end

"""

    masked_mean_lastdim(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the sliding-window mean of an array over the last dimension by using a mask.
If X is a CuArray, the computation takes place on the GPU.
"""
function masked_mean_lastdim(
    X::A,
    M::A,
) where {A<:Union{Matrix{T},CuArray{T,2}}} where {T<:Real}
    return (X * M) ./ reduce(+, M[:, 1])
end

"""

    var(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the variance of an array over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""
function var(X::A) where {A<:Union{Array{T},CuArray{T}}} where {T<:Real}
    Xmean = mean_lastdim(X)
    return var(X, Xmean)
end

function var(X::A, Xmean::A) where {A<:Union{Array{T},CuArray{T}}} where {T<:Real}
    lastdim = length(size(X))
    return reduce(+, (X .- Xmean) .^ 2, dims = lastdim) ./ T(size(X, lastdim) - 1)
end

"""

    masked_meansquare(X::A, M::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the mean squared value of a time series.
If the time series is detrended, this corresponds loosely to the variance.
If X is a CuArray, the computation takes place on the GPU.
"""
function masked_meansquare(
    X::A,
    M::A,
) where {A<:Union{Matrix{T},CuArray{T,2}}} where {T<:Real}
    return (X .^ 2) * M ./ reduce(+, M[:, 1])
end

"""

    skw(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the skewness of an array of dimension < 3 over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""
function skw(
    X::A,
    Xmean::A,
    Xvar::A,
) where {A<:Union{Array{T},CuArray{T}}} where {T<:AbstractFloat}
    lastdim = length(size(X))
    n = size(X, lastdim)
    cor = n^2 / (n - 1) / (n - 2)
    return cor .* (reduce(+, (X .- Xmean) .^ 3, dims = lastdim) ./ n ./ (Xvar .^ T(1.5)))
end

function skw(X::A) where {A<:Union{Array{T},CuArray{T}}} where {T<:AbstractFloat}
    Xmean = mean_lastdim(X)
    Xbvar = var(X)
    return skw(X, Xmean, Xbvar)
end

"""

    krt(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the kurtosis of an array of dimension < 3 over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""
function krt(
    X::A,
    Xmean::A,
    Xvar::A,
) where {A<:Union{Array{T},CuArray{T}}} where {T<:AbstractFloat}
    lastdim = length(size(X))
    n = size(X, lastdim)
    cor1 = (n + 1) * n / (n - 1) / (n - 2) / (n - 3)
    cor2 = 3 * (n - 1)^2 / (n - 2) / (n - 3)
    return cor1 .* (reduce(+, (X .- Xmean) .^ 4, dims = lastdim) ./ ((Xvar .^ T(2)))) .-
           cor2
end

function krt(X::A) where {A<:Union{Array{T},CuArray{T}}} where {T<:AbstractFloat}
    Xmean = mean_lastdim(X)
    Xvar = var(X)
    return krt(X, Xmean, Xvar)
end

#####################################################
# Analytic regression models
#####################################################
"""

    ar1_whitenoise(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the AR1 coefficient of an array of dimension < 3 over the last dimension.
Relies on the analytic solution of the least-square parameter estimation which reads:

If X is a CuArray, the computation takes place on the GPU.
"""
function ar1_whitenoise(x::Vector{T}) where {T<:Real}
    return (x[2:end]' * x[1:end-1]) / (x[1:end-1]' * x[1:end-1])
end

# AR1 coefficients for a vertical stack of residuals X and white noise assumption.
function ar1_whitenoise(X::A) where {A<:Union{Matrix{T},CuArray{T,2}}} where {T<:Real}
    lastdim = length(size(X))
    return reduce(+, X[:, 2:end] .* X[:, 1:end-1], dims = 2) ./
           reduce(+, X[:, 1:end-1] .* X[:, 1:end-1], dims = 2)
end

"""

    masked_ar1_whitenoise(X::A, M::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the same as [`ar1_ar1noise`](@ref ar1_whitenoise) but by means of a masking matrix.
This provides a significant speed-up for large data sets that are computed on GPU.
"""
function masked_ar1_whitenoise(
    X::A,
    M::A,
) where {A<:Union{SparseMatrixCSC{T,Int},CuArray{T,2}}} where {T<:Real}
    return ((X[:, 2:end] .* X[:, 1:end-1]) * M[1:end-1, :]) ./
           ((X[:, 1:end-1] .* X[:, 1:end-1]) * M[1:end-1, :])
end

# TODO: implement local Hurst exponent

#####################################################
# Analytic regression with correlated noise
#####################################################

function ar1_ar1noise() end

#####################################################
# Analytic regression uneven time spacing
#####################################################

function ar1_uneven_tspacing() end

#####################################################
# Numerical regression
#####################################################

function restoring_rate() end

function restoring_rate_gls() end

#####################################################
# Frequency spectrum
#####################################################

"""

    lfps(X::A; q_lowfreq=0.1) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the low-frequency power spectrum of a matrix over the 2nd dimension.
If X is a CuArray, the computation takes place on the GPU.
"""
function lfps(X::Matrix{T}; q_lowfreq = 0.1::AbstractFloat) where {T<:Real}
    P = abs.(rfft(X, 2))
    Pnorm = P ./ reduce(+, P, dims = 2)
    return reduce(+, Pnorm[:, 1:roundint(q_lowfreq * size(Pnorm, 2))], dims = 2)
end

function lfps(X::CuArray{T,2}; q_lowfreq = 0.1::AbstractFloat) where {T<:Real}
    P = abs.(CUDA.CUFFT.rfft(X, 2))
    Pnorm = P ./ reduce(+, P, dims = 2)
    return reduce(+, Pnorm[:, 1:roundint(q_lowfreq * size(Pnorm, 2))], dims = 2)
end
# Here multiple dispatch is needed as FFT functions are different

# TODO: implement wavelet analysis?

#####################################################
# Spatial indicators
#####################################################

# check out https://www.early-warning-signals.org/
function network_connectivity() end

function spatial_variance() end

function spatial_ar1() end

function recovery_length() end

function density_ratio() end
