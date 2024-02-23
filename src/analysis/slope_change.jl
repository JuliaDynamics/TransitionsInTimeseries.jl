import LsqFit

"""
    SlopeChangeConfig <: ChangesConfig
    SlopeChangeConfig(; indicator = nothing, kw...)

A configuration that can be given to [`estimate_changes`](@ref).
It estimates a change of slope in the timeseries by fitting two
connected linear segments to the timeseries,
returning the results (i.e., the two-linear fits) as [`SlopeChangeResults`](@ref).

## Keyword arguments
- indicator = nothing: if not nothing. Otherwise it should be a function f(x) -> Real.
  The slope fitting is then done over an indicator of the timeseries, which itself
  is estimated via a sliding window exactly as in [`SlidingWindowConfig`](@ref).
- `width_ind, stride_ind, whichtime`: exactly as in [`SlidingWindowConfig`](@ref)
  if `indicator` is not `nothing`.
"""
@kwdef struct SlopeChangeConfig{I, W} <: ChangesConfig
    indicators::I = nothing
    width_ind::Int = 100
    stride_ind::Int = 1
    whichtime::W = midpoint
end

function estimate_changes(config::ChangesConfig, x, t = eachindex(x))
    indicators = config.indicators
    # initialize time vectors
    if isnothing(indicators)
        # Skip indicators if they are nothing
        t_indicator = t
        x_indicator = x
    else
        t_indicator = windowmap(config.whichtime, t;
            width = config.width_ind, stride = config.stride_ind
        )
        x_indicator = windowmap(config.indicators, x;
            width = config.width_ind, stride = config.stride_ind
        )
    end
    p0 = guess_initial_p(x, t)
    fit = LsqFit.curve_fit(twolinear, t_indicator, x_indicator, p0)
    pbest = LsqFit.coef(fit)
    a, b, c, d = pbest
    t_change = (c - a)/(b - d)
    return SlopeChangeResults(t, x, t_indicator, x_indicator, [t_change], pbest, config, fit)
end

function guess_initial_p(x, t)
    midindex = (firstindex(x) + lastindex(x))÷2
    x1 = x[firstindex(x):midindex]
    x2 = x[midindex:lastindex(x)]
    t1 = t[firstindex(x):midindex]
    t2 = t[midindex:lastindex(x)]
    a, b = linreg(t1, x1)
    c, d = linreg(t2, x2)
    return [a, b, c, d]
end

import Statistics
function linreg(x, y)
    mx = Statistics.mean(x)
    my = Statistics.mean(y)
    b = Statistics.covm(x, mx, y, my)/Statistics.varm(x, mx)
    a = my - b*mx
    return a, b
end

function twolinear(t, p)
    a, b, c, d = p
    tcrit = (c - a)/(b - d)
    return @. ifelse(t < tcrit, a + b*t, c + d*t)
end

"""
    SlopeChangeResults <: ChangesResults

A struct containing the output of [`estimate_changes`](@ref) used with
[`SlopeChangeConfig`](@ref). It can be used for further analysis, visualization,
or given to [`significant_transitions`](@ref). The only significance type
that you can use this with [`significant_transitions`](@ref) is
[`SlopeChangeSignificance`](@ref).

It has the following fields that the user may access:

- `x`: the input timeseries.
- `t`: the time vector of the input timeseries.
- `x_indicator`, the indicator timeseries.
- `t_indicator`, the time vector of the indicator timeseries.
- `t_change`, the time the slope changes.
- `fitparams = a, b, c, d`, the fitted linear coefficients, `a + b*t` before.
  `t_change` and `c + d*t` after `t_change`.
"""
struct SlopeChangeResults{T, X, W, L} <: ChangesResults
    t # we don't parameterize these; they are only used for plotting
    x
    t_indicator::T
    x_indicator::X
    t_change::Vector{Float64}
    fitparams::Vector{Float64}
    config::W
    lsqfit::L
end
