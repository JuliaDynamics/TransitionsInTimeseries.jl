struct IndicatorEvolutionResults{T<:AbstractFloat}
    t_indicator::AbstractVector{T}
    x_indicator::AbstractVector{T}
    t_evolution::AbstractVector{T}
    x_evolution::AbstractVector{T}
    surr_x_evolution::AbstractMatrix{T}
end

#####################################################
# Trend metrics
#####################################################

@inline function compute_trend(
    t::Vector{T},
    x::Vector{T},
    indicator::Function,
    n_surr::Int,
    surr_method::S,
    trend_metric::Function,
    wv_indicator_width::Int,
    wv_indicator_stride::Int,
    wv_trend_width::Int,
    wv_trend_stride::Int,
) where {T<:AbstractFloat, S<:Surrogate}

    wv_indicator = WindowViewer(x, wv_indicator_width, wv_indicator_stride)
    t_indicator = t[wv_indicator.strided_indices]
    x_indicator = map(indicator, wv_indicator)

    wv_trend = WindowViewer(x_indicator, wv_trend_width, wv_trend_stride)
    t_evolution = t_indicator[wv_trend.strided_indices]
    x_evolution = map(trend_metric, wv_trend)

    sgen = surrogenerator(x, surr_method)
    surr_evolutions = fill(T(0.0), n_surr, length(x_evolution))
    for i in 1:n_surr
        s = sgen()
        wv_surr_indicator = WindowViewer(s, wv_indicator_width, wv_indicator_stride)
        surr_x_indicator = map(indicator, wv_surr_indicator)
        wv_surr_trend = WindowViewer(surr_x_indicator, wv_trend_width, wv_trend_stride)
        surr_trend_ts = map(trend_metric, wv_surr_trend)
        surr_evolutions[i, :] = surr_trend_ts
    end

    return IndicatorEvolutionResults(
        t_indicator,
        x_indicator,
        t_evolution,
        x_evolution,
        surr_evolutions,
    )
end

"""

    ridge_regression(x, y)

Perform ridge regression of `y` over `x` with regularization parameter `lambda`.
Return vector containing slope and offset.
If `lambda = 0`, linear regression is recovered (default case).
For more information, visit: https://en.wikipedia.org/wiki/Ridge_regression
"""
function ridge_regression(
    x::AbstractVector{T},
    y::AbstractVector{T};
    lambda::T = T(0),
) where {T<:Real}
    X = hcat(x, ones(length(x)))'
    eye = diagm(fill(T(1), 2))
    return inv(X * X' + lambda .* eye ) * X * y
end


"""

    ridge_regression_slope(x, y)

Extract the slope of ridge regression of `y` over `x`.
"""
function ridge_regression_slope(
    x::AbstractVector{T},
    y::AbstractVector{T};
    lambda::T = T(0),
) where {T<:Real}
    return ridge_regression(x, y, lambda = lambda)[1]
end

function ridge_regression_slope(
    y::AbstractVector{T}
) where {T<:Real}
    return ridge_regression(T.(eachindex(y)), y)[1]
end
#=> 
Trend metrics reexported by StatsBase:
corspearman(x, y=x)
corkendall(x, y=x)
<=#

#####################################################
# Change-point metrics
#####################################################