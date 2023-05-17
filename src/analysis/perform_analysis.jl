"""
    indicators_analysis(t, x, indicators::IndicatorsConfig, sigconfig::SignificanceConfig)

Perform an analysis of transition indicators for input timeseries `x` with time vector `t`
by specifying which indicators to use and with what sliding window ([`IndicatorsConfig`](@ref)),
as well as how to measure significant changes in the indicators ([`SignificanceConfig`](@ref)).
If `t` is not provided, it is assumed that `t=eachindex(x).`

Return the output as [`IndicatorsResults`](@ref).

This function performs the analysis described in the documentation
Example sections. It computes various indicators over sliding windows of `x`,
and then computes change metrics of over sliding windows of the indicators.
It does exactly the same computations for surrogates of `x` and use this
result to check for significance of the change metrics by computing the p-value.
The returned output contains all these computed timeseries.
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
    pval = zeros(X, len_change, n_ind)
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
        @inbounds for k in 1:sigconfig.n_surrogates
            s = sgen()
            windowmap!(indconfig.indicators[i], indicator_dummy, s;
                width = indconfig.width, stride = indconfig.stride)
            windowmap!(sigconfig.change_metrics[i_metric], change_dummy, indicator_dummy;
                width = sigconfig.width, stride = sigconfig.stride
            )
            # This should be replaced by a simple call of pvalue() in future.
            # However, the use of pvalue() is less trivial for an incremental
            # computation over the surrogates.
            if sigconfig.tail == :right
                pval[:, i] += c .< change_dummy
            elseif sigconfig.tail == :left
                pval[:, i] += c .> change_dummy
            else
                pr = c .< change_dummy
                pl = c .> change_dummy
                pval[:, i] += 2min.(pr, pl)
            end
        end
    end
    # pvalue = count / n_surrogates.
    pval ./= sigconfig.n_surrogates

    # put everything together in the output type
    return IndicatorsResults(
        t, x,
        indconfig.indicators, t_indicator, x_indicator,
        sigconfig.change_metrics, t_change, x_change, pval,
        sigconfig.surrogate_method,
    )
end