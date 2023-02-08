"""

    measure_significance(res, significance_metric)

Compute the `significance_metric` on the results `res` of an indicator evolution
analysis.
"""
@inline function measure_significance(
    res::IndicatorEvolutionResults{T},
    significance_metric::Function,
) where{T<:AbstractFloat}
    significances = similar(res.x_evolution)
    for j in eachindex(significances)
        significances[j] = significance_metric(x[j], res.surr_evolutions[:, j])
    end
end

"""

    get_quantile_idx(s)

Compute the indices of the quantile values within the vector `s`.
"""
@inline function get_quantile_idx(
    s::AbstractVector{T};
    q::T = T(0.95),
    symmetric::Bool = false,
) where{T<:AbstractFloat}
    sorted_idx = sortperm(s)
    i_lo = symmetric ? sorted_idx[roundint(n * (1 - q))] : 1
    i_hi = sorted_idx[roundint(n * q)]
    return i_lo, i_hi
end

"""

    get_normalized_interquantile(x, s)

Compute the normalized distance between an input value `x` and the
quantiles of a vector `s`.
"""
@inline function get_normalized_interquantile(
    x::T,
    s::AbstractVector{T};
    q::T = T(0.95),
    symmetric::Bool = false,
) where{T<:AbstractFloat}
    i_lo, i_hi = get_quantile_idx(s, q=q, symmetric=symmetric)
    s_mi = 0.5 * (s[i_lo] + s[i_hi])
    return (x - s_mi) / s_mi
end

roundint(x::Real) = round(Int(x))