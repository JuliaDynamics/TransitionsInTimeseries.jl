"""
QuantileSignificance(; p = 0.95, tail = :right) <: TransitionsSignificance

A configuration struct for significance testing [`significant_transitions`](@ref).
When used with [`WindowedIndicatorResults`](@ref), significance is estimated as
by comparing the value of each change metric with its `p`-quantile.
Values that exceed the `p`-quantile (if `tail = :right`)
or subseed the `1-p`-quantile (if `tail = :left`)
are deemed significant.
If `tail = :both` then either condition is checked.
"""
Base.@kwdef struct QuantileSignificance{P<:Real}
    p::P = 0.95
    tail::Symbol = :right
end

using Statistics: quantile

function significant_transitions(res::WindowedIndicatorResults, signif::QuantileSignificance)
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

