#####################################################
# Statistical moments
#####################################################

# TODO thorough tests for each function.

"""
    mean(X::AbstractArray)

Computes the mean of an array of dimension < 3 over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""

function mean(x::Vector{T}) where {T<:Real}
    return StatsBase.mean(x)
end

function mean(x::Matrix{T}) where {T<:Real}
    return StatsBase.mean(x, dims=2)
end

function mean(X::CuArray{T, 2}) where {T<:Real}
    return reduce( +, X, dims=2) ./ T(size(X, 2))
end

"""
    var(X::AbstractArray)

Computes the variance of an array of dimension < 3 over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""

function var(x::Vector{T}) where {T<:Real}
    return StatsBase.var(x)
end

function var(x::Matrix{T}) where {T<:Real}
    return StatsBase.var(x, dims=2)
end

# Unbiased estimator (consistent with StatsBase).
function var(X::CuArray{T, 2}, x_mean::CuArray{T, 2}) where {T<:Real}
    return reduce( +, (X .- x_mean).^2, dims=2) ./ T(size(X, 2) - 1)
end

function var(X::CuArray{T, 2}) where {T<:Real}
    x_mean = mean(X)
    return var(X, x_mean)
end

"""
    skw(X::AbstractArray)

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

"""
    krt(X::AbstractArray)

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

#####################################################
# Regression models
#####################################################

"""
    ar1_whitenoise(X::AbstractArray)

Computes the AR1 coefficient of an array of dimension < 3 over the last dimension.
If X is a CuArray, the computation takes place on the GPU.
"""

# AR1 coefficient of a vector x for white noise assumption.
# M. Mudelsee, Climate Time Series Analysis, eq 2.4
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

# On example nt = 1000, ns = 10000: speedup is factor 10 compared to GPU + CPU loop.
function ar1_whitenoise_acc(X::CuArray{T, 2}, M::CuArray{T, 2}) where {T<:Real}
    return ((X[:, 2:end] .* X[:, 1:end-1]) * M[1:end-1, :]) ./ ((X[:, 1:end-1] .* X[:, 1:end-1]) * M[1:end-1, :])
end

# TODO implement the TIs below.
function ar1_ar1noise()
end

function ar1_uneven_tspacing()
end

function restoring_rate()
end

function restoring_rate_gls()
end

#####################################################
# Frequency spectrum
#####################################################

"""
    lfps(X::AbstractArray; q_lowfreq=0.1)

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

#####################################################
# Spatial identifiers
#####################################################

# TODO implement EWSs below
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