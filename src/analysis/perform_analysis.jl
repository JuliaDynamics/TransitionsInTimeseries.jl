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
    if sigconfig.tail == :both
        pval_right = zeros(X, len_change, n_ind)
        pval_left = zeros(X, len_change, n_ind)
    end
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
            elseif sigconfig.tail == :both
                pval_right[:, i] += c .< change_dummy
                pval_left[:, i] += c .> change_dummy
            end
        end
        if sigconfig.tail == :both
            pval[:, i] .= 2min.(pval_right[:, i], pval_left[:, i])
        end
    end
    pval ./= sigconfig.n_surrogates

    # put everything together in the output type
    return IndicatorsResults(
        t, x,
        indconfig.indicators, t_indicator, x_indicator,
        sigconfig.change_metrics, t_change, x_change, pval,
        sigconfig.surrogate_method,
    )
end

"""

    transition_flags(results::IndicatorsResults, p_threshold::Real)

Return `tflags_indicators::Vector{Vector}` and `tflags_andicators::Vector`.
The former contains the time steps at which the p-value of an indicator change is
below `p_threshold` (further called flags). The flags of the `i`-th indicator can
be obtained by calling tflags_indicators[i]. The latter contains the time steps
when all p-values of all computed indicator changes are synchronously below `p_threshold`.
"""
function transition_flags(results::IndicatorsResults, p_threshold::Real)
    if p_threshold >= 1 || p_threshold <= 0
        error("Threshold must be a value between 0 and 1.")
    end

    thresholded_p = results.pval .< threshold
    tflags_indicators = [t[thresholded_p[:, j]] for j in axes(thresholded_p, 2)]
    tflags_andicators = results.t_change[reduce(&, thresholded_p, dims=2)]
    return tflags_indicators, tflags_andicators
end