"""

    measure_significance(res, significance_metric)

Compute the `significance_metric` on the results `res` of an indicator evolution
analysis.
"""
@inline function measure_significance(
    res::IndicatorEvolutionResults{T},
    significance_metric::Function,
) where {T<:AbstractFloat}
    return measure_significance(
        res.x_evolution,
        res.surr_x_evolution,
        significance_metric,
    )
end

@inline function measure_significance(
    x_evolution::AbstractVector{T},
    surr_x_evolution::AbstractMatrix{T},
    significance_metric::Function,
) where {T<:AbstractFloat}
    significance = similar(x_evolution)
    for j in eachindex(significance)
        significance[j] = significance_metric(
            x_evolution[j],
            view(surr_x_evolution, :, j),
        )
    end
    return significance
end

"""

    get_quantile_idx(s)

Compute the indices of the quantile values within the vector `s`.
"""
@inline function quantile_idx(
    s::AbstractVector{T};
    q::T = T(0.95),
    symmetric::Bool = false,
) where{T<:AbstractFloat}
    n = length(s)
    sorted_idx = sortperm(s)
    i_lo = symmetric ? sorted_idx[intround(n * (1 - q))] : sorted_idx[1]
    i_hi = sorted_idx[intround(n * q)]
    return i_lo, i_hi
end

"""

    normalized_quantile_distance(x, s)

Compute the normalized distance between an input value `x` and the
quantiles of a vector `s`. 
"""
@inline function normalized_quantile_distance(
    x::T,
    s::AbstractVector{T};
    q::T = T(0.95),
    symmetric::Bool = false,
) where{T<:AbstractFloat}
    i_lo, i_hi = quantile_idx(s, q=q, symmetric=symmetric)
    center = 0.5 * (s[i_lo] + s[i_hi])
    interq = 0.5 * (s[i_hi] - s[i_lo])
    return (x - center) / abs(interq)
end

"""

    gaussian_quantile(x, s)

Compute the normalized distance between an input value `x` and the
quantiles of a vector `s`. 
"""
@inline function gaussian_quantile(
    x::T,
    s::AbstractVector{T};
    nstd::T=2.0,
) where{T<:AbstractFloat}
    return T( (x - mean(s)) / (nstd * std(s)) )
end

"""

    which_quantile(x, s)

Compute the quantile represented by the input value `x` w.r.t. the vector `s`. 
"""
@inline function which_quantile(
    x::T,
    s::AbstractVector{T},
) where{T<:AbstractFloat}
    return sum(x .> s) / length(s)
end

intround(x::Real) = Int(round(x))