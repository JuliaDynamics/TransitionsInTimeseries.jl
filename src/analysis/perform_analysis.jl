"""
    indicators_analysis(x, indicators::IndicatorsConfig, sigconfig::SignificanceConfig)

Perform an analysis of transition indicators for input timeseries `x`
by specifying which indicators to use and with what sliding window
([`IndicatorsConfig`](@ref)), as well as how to measure significant changes in the
indicators ([`SignificanceConfig`](@redf)).

Return the output as [`IndicatorsResults`](@ref).

This function performs the analysis described in the documentation
Example sections. It computes various indicators over sliding windows of `x`,
and then computes change metrics of over sliding windows of the indicators.
In parallel it does exactly the same computations for surrogates of `x`.
The returned output contains all these computed timeseries and can be given
to [`indicators_significance`](@ref).
"""
function indicators_analysis(x, indconfig::IndicatorsConfig, sigconfig::SignificanceConfig)
    t = eachindex(x)
    return indicators_analysis(t, x, indconfig, sigconfig)
end

function indicators_analysis(t::AbstractVector, x, indconfig::IndicatorsConfig, sigconfig::SignificanceConfig)
    X = eltype(x)
    # initialize sizes
    n_ind = length(indconfig.indicators)
    if length(sigconfig.change_metrics) == n_ind
        one2one = true
    elseif length(sigconfig.change_metrics) == 1
        one2one = false
    else
        error("The amount of change metrics must either be 1 or be the same " *
        "as the amount of indicators.")
    end
    t_indicator = indconfig.t_indicator
    t_change = sigconfig.t_change
    len_ind = length(t_indicator)
    len_change = length(t_change)
    # initialize array containers
    x_indicator = zeros(X, len_ind, n_ind)
    x_change = zeros(X, len_change, n_ind)
    s_change = zeros(X, len_change, n_ind, sigconfig.n_surrogates)
    indicator_dummy = zeros(X, len_ind)
    change_dummy = zeros(X, len_change)
    sgen = surrogenerator(x, sigconfig.surrogate_method, sigconfig.rng)
    # Actual computations
    @inbounds for i in 1:n_ind
        i_metric = one2one ? i : 1
        # indicator timeseries
        z = view(x_indicator, :, i)
        windowmap!(indconfig.indicators[i], z, x;
            width = indconfig.width, stride = indconfig.stride)
        # change metric timeseries
        c = view(x_change, :, i)
        windowmap!(sigconfig.change_metrics[i_metric], c, z;
            width = sigconfig.width, stride = sigconfig.stride)
        # surrogates
        for k in 1:sigconfig.n_surrogates
            s = sgen()
            windowmap!(indconfig.indicators[i], indicator_dummy, s;
                width = indconfig.width, stride = indconfig.stride)
            windowmap!(sigconfig.change_metrics[i_metric], change_dummy, indicator_dummy;
                width = sigconfig.width, stride = sigconfig.stride
            )
            s_change[:, i, k] .= change_dummy
        end
    end
    # put everything together in the output type
    return IndicatorsResults(
        t, x,
        indconfig.indicators, t_indicator, x_indicator,
        sigconfig.change_metrics, t_change, x_change, s_change,
        sigconfig.surrogate_method,
    )
end