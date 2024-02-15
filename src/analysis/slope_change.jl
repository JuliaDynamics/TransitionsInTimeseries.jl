import LsqFit

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
    return SlopeChangeResults(t, x, t_indicator, x_indicator, t_change, pbest, config, fit)
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

struct SlopeChangeResults{T, X, W, L} <: ChangesResults
    t # we don't parameterize these; they are only used for plotting
    x
    t_indicator::T
    x_indicator::X
    t_change::Float64
    fitparams::Vector{Float64}
    config::W
    lsqfit::L
end
