#####################################################
# ::Params and <:Params
#####################################################

abstract type Params end
abstract type StructFunction <: Function end

"""
    IndicatorsParams <: Params

A `struct` containing the parameters inherent to the functions estimating the indicators.
`IndicatorsParams` comes with default choices and can simply be initialized by
`IndicatorsParams()`. The parameters can be user-defined through the keyword arguments.

## Keyword arguments
- `q_lofreq::Real`: between `0` and `1`, defines the fraction of low frequencies
  in the power density spectrum.
- `m_perm::Int`: defines the order of symbolic permutations for the permutation entropy.
"""
Base.@kwdef struct IndicatorsParams <: Params
    q_lofreq::Real = 0.1
    m_perm::Int = 3
end

"""
    ChangeMetricsParams <: Params

A `struct` containing the parameters inherent to the functions estimating the change
metrics. `ChangeMetricsParams` comes with default choices and can simply be initialized
by `IndicatorsParams()`. Parameters can be user-defined through the keyword arguments.

## Keyword arguments
- `lambda_ridge::Real`: between `0` and `Inf`, the regularization parameter of the
  ridge regression.
"""
Base.@kwdef struct ChangeMetricsParams <: Params
    lambda_ridge::Real = 0.0
end

"""
    init_metrics(metrics, t, p)

Initialize the functions contained in `metrics::Vector` that require an initialization
based on the time vector `t:.AbstractVector` and parameters `p::Params`.
These functions are typed as `StructFunction` and share the same initialization
syntax `metric(t, p)`.

The output `f::Vector{Function}` contains the ready-to-use metrics.
"""
function init_metrics(
    metrics::Vector,
    t::AbstractVector,
    p::Params,
)
    f = Function[]
    for metric in metrics
        if metric in subtypes(StructFunction)
            push!(f, metric(t, p))
        elseif isa(metric, Function)
            push!(f, metric)
        else
            error("The function $metric has an erroneous type.")
        end
    end
    return f
end

#####################################################
# Default choices
#####################################################

default_n_surrogates() = 10_000
default_surrogate_method() = RandomFourier()

function default_window_width(x)
    w = max(length(x)÷100, 10)
    if w < 20
        @warn "The window width chosen by default is w = $w and might be too small!"
    else
        @info "The window width chosen by default is w = $w."
    end
    return w
end

function default_window_stride(x)
    s = 1
    @info "The window stride chosen by default is s = $s."
    return s
end