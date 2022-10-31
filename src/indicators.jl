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

"""@docs
    mean(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the mean of an array over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""
function mean(X::Array{T}) where {T<:Real}
    return StatsBase.mean(X, dims=length(size(X)))
end

function mean(X::CuArray{T, 2}) where {T<:Real}
    reduce_dim = length(size(X))
    return reduce( +, X, dims=reduce_dim) ./ T(size(X, reduce_dim))
end

"""@docs
    masked_mean(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the sliding-window mean of an array over the last dimension by using a mask.
If X is a CuArray, the computation takes place on the GPU.
"""
function masked_mean(X::A, M::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}} where {T<:Real}
    return (X * M) ./ reduce(+, M[:, 1])
end

"""@docs
    var(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the variance of an array over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""
function var(X::Array{T}) where {T<:Real}
    return StatsBase.var(X, dims=length(size(X)))
end

# Unbiased estimator (consistent with StatsBase).
function var(X::CuArray{T, 2}, x_mean::CuArray{T, 2}) where {T<:Real}
    return reduce( +, (X .- x_mean).^2, dims=2) ./ T(size(X, 2) - 1)
end

function var(X::CuArray{T, 2}) where {T<:Real}
    x_mean = mean(X)
    return var(X, x_mean)
end

"""@docs
    masked_meansquare(X::A, M::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the mean squared value of a time series.
If the time series is detrended, this corresponds loosely to the variance.
If X is a CuArray, the computation takes place on the GPU.
"""
function mean_square(X::A, M::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}} where {T<:Real}
    return (X .^ 2) * M ./ reduce(+, M[:, 1])
end

"""@docs
    skw(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the skewness of an array of dimension < 3 over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""
function skw(x::Vector{T}) where {T<:Real}
    return StatsBase.skewness(x)
end

function skw(X::Matrix{T}, x_mean::Matrix{T}, x_var::Matrix{T}) where {T<:Real}
    return reduce( +, (X .- x_mean).^3, dims=2) ./ size(X, 2) ./ (x_var .^ T(1.5))
end

function skw(X::Matrix{T}) where {T<:Real}
    x_mean = StatsBase.mean(X, dims=2)
    x_var = var(X)
    return skw(X, x_mean, x_var)
end

function skw(X::CuArray{T, 2}, x_mean::CuArray{T, 2}, x_var::CuArray{T, 2}) where {T<:Real}
    return reduce( +, (X .- x_mean).^3, dims=2) ./ size(X, 2) ./ (x_var .^ T(1.5))
end

function skw(X::CuArray{T, 2}) where {T<:Real}
    x_mean = mean(X)
    x_var = var(X)
    return skw(X, x_mean, x_var)
end
# TODO compare to skw of StatsBase

"""@docs
    krt(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the kurtosis of an array of dimension < 3 over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""
function krt(x::Vector{T}) where {T<:Real}
    return StatsBase.kurtosis(x)
end

function krt(X::Matrix{T}, x_mean::Matrix{T}, x_var::Matrix{T}) where {T<:Real}
    return reduce( +, (X .- x_mean).^3, dims=2) ./ size(X, 2) ./ (x_var .^ T(1.5))
end

function krt(X::Matrix{T}) where {T<:Real}
    x_mean = StatsBase.mean(X, dims=2)
    x_var = var(X)
    return reduce( +, (X .- x_mean).^3, dims=2) ./ size(X, 2) ./ (x_var .^ T(1.5))
end

function krt(X::CuArray{T, 2}, x_mean::CuArray{T, 2}, x_var::CuArray{T, 2}) where {T<:Real}
    return reduce( +, (X .- x_mean) .^ 4, dims=2 ) ./ size(X, 2) ./ (x_var .^ 2)
end

function krt(X::CuArray{T, 2}) where {T<:Real}
    x_mean = mean(X)
    x_var = var(X)
    return krt(X, x_mean, x_var)
end
# TODO compare to krt of StatsBase

#####################################################
# Analytic regression models
#####################################################
"""@docs
    ar1_whitenoise(X::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the AR1 coefficient of an array of dimension < 3 over the last dimension.
Relies on the analytic solution of the least-square parameter estimation which reads:

If X is a CuArray, the computation takes place on the GPU.
"""
function ar1_whitenoise(x::Vector{T}) where {T<:Real}
    return (x[2:end]' * x[1:end-1]) / (x[1:end-1]' * x[1:end-1])
end

# AR1 coefficients for a vertical stack of residuals X and white noise assumption.
function ar1_whitenoise(X::Array{T}) where {T<:Real}
    return reduce( +, X[:, 2:end] .* X[:, 1:end-1], dims=2) ./ reduce( +, X[:, 1:end-1] .* X[:, 1:end-1], dims=2)
end

# GPU accelerated AR1 coefficients for a vertical stack of time series X for white noise assumption.
function ar1_whitenoise(X::CuArray{T, 2}) where {T<:Real}
    return reduce( +, X[:, 2:end] .* X[:, 1:end-1], dims=2) ./ reduce( +, X[:, 1:end-1] .* X[:, 1:end-1], dims=2)
end

"""@docs
    masked_ar1_whitenoise(X::A, M::A) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the same as [`ar1_ar1noise`](@ref ar1_whitenoise) but by means of a masking matrix.
This provides a significant speed-up for large data sets that are computed on GPU.
"""
function masked_ar1_whitenoise(X::CuArray{T, 2}, M::CuArray{T, 2}) where {T<:Real}
    return ((X[:, 2:end] .* X[:, 1:end-1]) * M[1:end-1, :]) ./ ((X[:, 1:end-1] .* X[:, 1:end-1]) * M[1:end-1, :])
end

# TODO: implement local Hurst exponent

#####################################################
# Analytic regression with correlated noise
#####################################################

function ar1_ar1noise()
end

#####################################################
# Analytic regression uneven time spacing
#####################################################

function ar1_uneven_tspacing()
end

#####################################################
# Numerical regression
#####################################################

function restoring_rate()
end

function restoring_rate_gls()
end

#####################################################
# Frequency spectrum
#####################################################

"""@docs
    lfps(X::A; q_lowfreq=0.1) where {A<:Union{Matrix{T}, CuArray{T, 2}}}

Computes the low-frequency power spectrum of a matrix over the 2nd dimension.
If X is a CuArray, the computation takes place on the GPU.
"""
function lfps(X::Matrix{T}; q_lowfreq=0.1::AbstractFloat) where {T<:Real}
    P = abs.(rfft(X, 2))
    Pnorm = P ./ reduce(+, P, dims=2)
    return reduce(+, Pnorm[:, 1:roundint(q_lowfreq * size(Pnorm, 2))], dims=2)
end

function lfps(X::CuArray{T, 2}; q_lowfreq=0.1::AbstractFloat) where {T<:Real}
    P = abs.(CUDA.CUFFT.rfft( X, 2 ))
    Pnorm = P ./ reduce(+, P, dims=2)
    return reduce(+, Pnorm[:, 1:roundint(q_lowfreq * size(Pnorm, 2))], dims=2)
end

# TODO: implement wavelet coefficient

#####################################################
# Spatial indicators
#####################################################

# check out https://www.early-warning-signals.org/
function network_connectivity()
end

function spatial_variance()
end

function spatial_ar1()
end

function recovery_length()
end

function density_ratio()
end