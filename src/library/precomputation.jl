abstract type PrecomputableFunction end

"""
    precompute_metrics(metrics, t)

Precompute functions contained in `metrics::Vector` that can be precomputed
based on the time vector `t::AbstractVector`. These functions are typed as
`PrecomputableFunction` and share the same initialization syntax
`precompute(metric, t)`.

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

