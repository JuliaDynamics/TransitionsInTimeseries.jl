"""
    indicators_analysis(x, indicators::IndicatorsConfig, significance::SignificanceConfig)

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
function indicators_analysis(x, indicators::IndicatorsConfig, significance::SignificanceConfig)
    # TODO: Allow `x` to be a `Timeseries`, that has its own time vector.
    t = eachindex(x)
    X = eltype(x)
    # initialize sizes
    n_ind = length(indicators.indicators)
    if length(significance.change_metrics) == n_ind
        one2one = true
    elseif length(significance.change_metrics) == 1
        one2one = false
    else
        error("The amount of change metrics must either be 1 or be the same " *
        "as the amount of indicators.")
    end
    len_ind = length(WindowViewer(x; indicators.window_kwargs...))
    len_change = length(WindowViewer(1:len_ind; width = significance.width, stride = significance.stride))
    # initialize array containers
    x_indicator = zeros(X, len_ind, n_ind)
    x_change = zeros(X, len_change, n_ind)
    s_change = zeros(X, len_change, n_ind, significance.n_surrogates)
    t_indicator = windowmap(midpoint, t; indicators.window_kwargs...)
    t_change = windowmap(midpoint, t_indicator; width = significance.width, stride = significance.stride)
    indicator_dummy = zeros(X, len_ind)
    change_dummy = zeros(X, len_change)
    sgen = surrogenerator(x, significance.surrogate_method, significance.rng)
    # Actual computations
    @inbounds for i in 1:n_ind
        i_metric = one2one ? i : 1
        # indicator timeseries
        z = view(x_indicator, :, i)
        windowmap!(indicators.indicators[i], z, x; indicators.window_kwargs...)
        # change metric timeseries
        c = view(x_change, :, i)
        windowmap!(significance.change_metrics[i_metric], c, z;
            width = significance.width, stride = significance.stride
        )
        # surrogates
        for k in 1:significance.n_surrogates
            s = sgen()
            windowmap!(indicators.indicators[i], indicator_dummy, s; indicators.window_kwargs...)
            windowmap!(significance.change_metrics[i_metric], change_dummy, indicator_dummy;
                width = significance.width, stride = significance.stride
            )
            s_change[:, i, k] .= change_dummy
        end
    end
    # put everything together in the output type
    return IndicatorsResults(
        t, x,
        indicators.indicators, t_indicator, x_indicator,
        significance.change_metrics, t_change, x_change, s_change,
        significance.surrogate_method,
    )
end
