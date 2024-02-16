"""
    ThresholdSignificance(threshold::Real; tail = :right) <: TransitionsSignificance

A configuration struct for significance testing in [`significant_transitions`](@ref).
Significance is estimated by comparing the value of each change metric with
the given `threshold`.
Values that exceed the threshold (if `tail = :right`)
or subseed the threshold (if `tail = :left`)
are deemed significant.
If `tail = :both` then either condition is checked.
"""
struct ThresholdSignificance{T<:Real}
    # TODO: threshold<:Real does not make sense. We should have as many thresholds as metrics.
    # TODO: once fixed, change the tests accordingly
    threshold::T
    tail::Symbol
end
ThresholdSignificance(threshold; tail = :right) = ThresholdSignificance(threshold, tail)

function significant_transitions(res::IndicatorsChangesResults, signif::ThresholdSignificance)
    flags = similar(res.x_change, Bool)
    t = signif.threshold
    for (i, x) in enumerate(eachcol(res.x_change))
        flag = view(flags, :, i)
        if signif.tail == :right
            @. flag = x > t
        elseif signif.tail == :left
            @. flag = x < t
        # FIXME: the case below does not make sense, since it always returns 1
        # It should rather be using two distinct thresholds
        elseif signif.tail == :both
            @. flag = (x < t) | (x > t)
        else
            error("`tail` can be only `:left, :right, :both`. Got $(tail).")
        end
    end
    return flags
end


"""
    QuantileSignificance(; p = 0.95, tail = :both) <: TransitionsSignificance

A configuration struct for significance testing [`significant_transitions`](@ref).
When used with [`IndicatorsChangesResults`](@ref), significance is estimated
by comparing the value of each change metric with its `p`-quantile.
Values that exceed the `p`-quantile (if `tail = :right`)
or subseed the `1-p`-quantile (if `tail = :left`)
are deemed significant.
If `tail = :both` then either condition is checked.

`QuantileSignficance` guarantees that some values will be significant by
the very definition of what a quantile is.
See also [`SigmaSignificance`](@ref) that is similar but does not have this guarantee.
"""
Base.@kwdef struct QuantileSignificance{P<:Real}
    p::P = 0.95
    tail::Symbol = :both
end

using Statistics: quantile

function significant_transitions(res::IndicatorsChangesResults, signif::QuantileSignificance)
    flags = similar(res.x_change, Bool)
    for (i, x) in enumerate(eachcol(res.x_change))
        qmin, qmax = quantile(x, (1 - signif.p, signif.p))
        flag = view(flags, :, i)
        if signif.tail == :right
            @. flag = x > qmax
        elseif signif.tail == :left
            @. flag = x < qmin
        elseif signif.tail == :both
            @. flag = (x < qmin) | (x > qmax)
        else
            error("`tail` can be only `:left, :right, :both`. Got $(tail).")
        end
    end
    return flags
end

"""
    SigmaSignificance(; factor = 3.0, tail = :both) <: TransitionsSignificance

A configuration struct for significance testing [`significant_transitions`](@ref).
When used with [`IndicatorsChangesResults`](@ref), significance is estimated
by comparing how many standard deviations (`σ`) the value exceeds the mean value (`μ`).
Values that exceed (if `tail = :right`) `μ + factor*σ`, or subseed (if `tail = :left`) `μ - factor*σ`
are deemed significant.
If `tail = :both` then either condition is checked.

`factor` can also be a vector of values,
in which case a different value is used for each change metric.

See also [`QuantileSignificance`](@ref).
"""
Base.@kwdef struct SigmaSignificance{P}
    factor::P = 0.95    # TODO: Default value should be an integer, Typically 1, 2 or 3.
    tail::Symbol = :both
end

using Statistics: std, mean

function significant_transitions(res::IndicatorsChangesResults, signif::SigmaSignificance)
    flags = similar(res.x_change, Bool)
    for (i, x) in enumerate(eachcol(res.x_change))
        μ = mean(x)
        σ = std(x; mean = μ)
        factor = signif.factor isa AbstractVector ? signif.factor[i] : signif.factor
        flag = view(flags, :, i)
        if signif.tail == :right
            @. flag = x > μ + factor*σ
        elseif signif.tail == :left
            @. flag = x < μ - factor*σ
        elseif signif.tail == :both
            @. flag = (x < μ - factor*σ) | (x > μ + factor*σ)
        else
            error("`tail` can be only `:left, :right, :both`. Got $(tail).")
        end
    end
    return flags
end

# TODO: compare2tail() should prevent the duplicated if-else checking the tails above