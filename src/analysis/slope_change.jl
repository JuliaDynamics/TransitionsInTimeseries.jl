import LsqFit

struct SlopeChangeConfig{I} <: ChangesConfig
    indicators::I
    width_ind::Int
    stride_ind::Int
    whichtime::W
end

struct SlopeChangeResults{T} <: ChangesResults
    t::T # original time vector
    x::X # original timeseries
    t_indicator::IT
    x_indicator::IX
    t_change::Float64 # only 1 change exists
    fitparams::Vector{Float64}
end

function estimate_changes(config::ChangesConfig, x, t = eachindex(x))
    indicators = config.indicators
    # initialize time vectors
    if isnothing(indicators)
        # Skip indicators if they are nothing
        t_indicator = t
    else
        t_indicator = windowmap(config.whichtime, t; width = config.width_ind,
            stride = config.stride_ind)
    end

    p0 = guess_initial_p(x, t)
    fit = LsqFit.curve_fit(twolinear, t_indicator, x_indicator, p0)
    pbest = LsqFit.coef(fit)

    t_change = (a, b, c, d = pbest; (c - a)/(b - d))
    return SlopeChangeResults(t, x, t_indicator, x_indicator, t_change, pbest)
end

function twolinear(t, p)
    a, b, c, d = p
    tcrit = (c - a)/(b - d)
    return @. ifelse(t < tcrit, a + b*t, c + d*t)
end
