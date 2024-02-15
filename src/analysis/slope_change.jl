import LsqFit

struct SlopeChangeConfig{I} <: ChangesConfig
    indicators::I
    width_ind::Int
    stride_ind::Int
    whichtime::W
end

struct SlopeChangeResults{
        TT, T<:Real, X<:Real, XX<:AbstractVector{X}, IT,
        W, L
    } <: ChangesResults
    t::T # original time vector
    x::X # original timeseries
    t_indicator::IT
    x_indicator::IX
    t_change::Float64
    fitparams::Vector{Float64}
    config::W
    lsqfit::L
end

function estimate_changes(config::ChangesConfig, x, t = eachindex(x))
    indicators = config.indicators
    # initialize time vectors
    if isnothing(indicators)
        # Skip indicators if they are nothing
        t_indicator = t
    else
        t_indicator = windowmap(config.whichtime, t;
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
