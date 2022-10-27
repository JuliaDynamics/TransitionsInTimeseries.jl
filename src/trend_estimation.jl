

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

# This implementation works fine for even time spacing
# because the time vector does not need to be shifted.
# TODO: make more general!

#####################################################
# Kendall tau
#####################################################

"""
    percentile_significance(ref_stat::AbstractArray, sur_stat::Matrix{T}, ns::Int, nx::Int)

Get percentile significance of reference (original time series) versus surrogate statistics.
"""
function kendall_tau(y::Vector{T}, t::Vector{T}) where {T<:Real}
    return StatsBase.corkendall(t, y)
end

function kendall_tau(Y::Matrix{T}, t::Vector{T}) where {T<:Real}
    return mapslices(x -> StatsBase.corkendall(t, x), Y, dims = 2)
end

function kendall_tau(Z::CuArray{T, 2}, t::Vector{T}) where {T<:Real}
    Y = Array(Z)
    return kendall_tau(Y, t)
end

# Scaling with mean (or other function) might be a better metric than
# just the trend! Typically: AR1 needs to be close to 1.
function scaled_kendall_tau(Y::Matrix{T}, t::Vector{T}; scaling::Function=StatsBase.mean) where {T<:Real}
    return mapslices(x -> StatsBase.corkendall(t, x) * scaling(x), Y, dims = 2)
end

# TODO: add Spearman’s ρ rank correlation or the Pearson’s correlation coefficient.

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

Slide trend estimation over computed indicator time series `X`. 
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
