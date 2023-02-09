struct IndicatorEvolutionResults{T<:AbstractFloat}
    t_indicator::AbstractVector{T}
    x_indicator::AbstractVector{T}
    t_evolution::AbstractVector{T}
    x_evolution::AbstractVector{T}
    surr_x_evolution::AbstractMatrix{T}
end

"""
    indicator_evolution(t, x, indicator, n_surr, surr_method, evolution_metric,
                        m, wv_indicator_width, wv_indicator_stride,
                        wv_evolution_width, wv_evolution_stride)


"""
@inline function indicator_evolution(
    t::AbstractVector{T},
    x::AbstractVector{T},
    indicator::Function,
    n_surr::Int,
    surr_method::S,
    evolution_metric::Function,
    wv_indicator_width::Int,
    wv_indicator_stride::Int,
    wv_evolution_width::Int,
    wv_evolution_stride::Int,
) where {T<:AbstractFloat, S<:Surrogate}

    wv_indicator = WindowViewer(x, wv_indicator_width, wv_indicator_stride)
    t_indicator = t[wv_indicator.strided_indices]
    x_indicator = map(indicator, wv_indicator)

    wv_evolution = WindowViewer(x_indicator, wv_evolution_width, wv_evolution_stride)
    t_evolution = t_indicator[wv_evolution.strided_indices]
    x_evolution = map(evolution_metric, wv_evolution)

    sgen = surrogenerator(x, surr_method)
    surr_evolutions = fill(T(0.0), n_surr, length(x_evolution))
    @inline for i in 1:n_surr
        s = sgen()
        wv_surr_indicator = WindowViewer(s, wv_indicator_width, wv_indicator_stride)
        surr_x_indicator = map(indicator, wv_surr_indicator)
        wv_surr_evolution = WindowViewer(surr_x_indicator, wv_evolution_width, wv_evolution_stride)
        surr_evolution_ts = map(evolution_metric, wv_surr_evolution)
        surr_evolutions[i, :] = surr_evolution_ts
    end

    return IndicatorEvolutionResults(
        t_indicator,
        x_indicator,
        t_evolution,
        x_evolution,
        surr_evolutions,
    )
end

#####################################################
# Trend metrics
#####################################################
"""

    ridge(t, x; lambda)

Perform ridge regression of `x` over `t` with regularization parameter `lambda`.
Return vector containing slope and offset.
If `lambda = 0`, linear regression is recovered (default case).
For more information, visit: https://en.wikipedia.org/wiki/Ridge_regression
"""
@inline function ridge(
    t::AbstractVector{T},
    x::AbstractVector{T};
    lambda::T = T(0),
) where {T<:Real}
    M = precompute_ridge(t, lambda = lambda)
    return ridge(M, x)
end

@inline function precompute_ridge(
    t::AbstractVector{T};
    lambda::T = T(0),
) where {T<:Real}
    TT = hcat(t, ones(T, length(t)))'
    return inv(TT * TT' + lambda .* I(2) ) * TT
end

@inline function ridge(
    M::AbstractMatrix{T},
    x::AbstractVector{T},
) where {T<:Real}
    return M * x
end

"""

    ridge_slope(x, y)

Extract the slope of ridge regression of `y` over `x`.
"""
function ridge_slope(
    x::AbstractVector{T},
    y::AbstractVector{T};
    lambda::T = T(0),
) where {T<:Real}
    return ridge(x, y, lambda = lambda)[1]
end

function ridge_slope(
    y::AbstractVector{T},
) where {T<:Real}
    return ridge(T.(eachindex(y)), y)[1]
end

"""

    precompute_ridge_slope(t)

For evenly spaced time series, precomputation can be performed to accelerate
the computation of the slope obtained by ridge regression over a time-span `t`.
"""
@inline function precompute_ridge_slope(
    t::AbstractVector{T};
    lambda::T = T(0),
) where {T<:Real}
    M = precompute_ridge(t, lambda=lambda)
    return M[1, :]
end

@inline function precomputed_ridge_slope(
    x::AbstractVector{T},
    m::AbstractVector{T},
) where {T<:AbstractFloat} 
    return m' * x
end



#=> 
Trend metrics reexported by StatsBase:
corspearman(x, y=x)
corkendall(x, y=x)
<=#

#####################################################
# Change-point metrics
#####################################################