"""

    PrecomputableFunction

Supertype of structs containing the necessary field to precompute a `::Function` by:

```julia
precompute(f::PrecomputableFunction, t)
```
"""
abstract type PrecomputableFunction end

"""
    precompute_metrics(metrics::Vector, t::AbstractVector)

Precompute functions contained in `metrics` that can be precomputed
based on the time vector `t`. Relies on abstract type `PrecomputableFunction`.

The output `f::Vector{Function}` contains the ready-to-use metrics.
"""
function precompute_metrics(metrics::Vector, t::AbstractVector)
    f = Function[]
    for metric in metrics
        if typeof(metric) in subtypes(PrecomputableFunction)
            push!(f, precompute(metric, t))
        elseif isa(metric, Function)
            push!(f, metric)
        else
            error("The function $metric has an erroneous type.")
        end
    end
    return f
end

