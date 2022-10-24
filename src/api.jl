"""
    compute_TIsignificance(
        t::Vector{T},
        X::CuArray{T, 2},
        pindctr::WindowingParams,
        pidtrend::WindowingParams,
        ns::Int,
        ti_function_list = [ar1_whitenoise, cuvar],
        trend_measure_function = ridge_regression_slope,
        window = centered_window;
        kwargs...,
    )

Compute various transition identifiers (provided in `ti_function_list`) for a residual `X` over time `t`.
"""

function TIs_significance(
    t::Vector{T},
    X::CuArray{T,2},
    pindctr::WindowingParams,
    pidtrend::WindowingParams,
    ns::Int,
    ti::Vector{Function},
    trend_measure::Function = ridge_regression_slope,
    window::Function = centered_window;
    kwargs...,
) where {T}

    nx, nt = size(X)
    S = generate_stacked_fourier_surrogates(X, ns)
    ti_labels = string.(ti_function_list)

    tindctr = trim_wndw(t, pindctr, window)
    tidtrend = trim_wndw(tindctr, pidtrend, window)

    reference_ti = Dict{String,CuArray{T,2}}()
    surrogate_ti = Dict{String,CuArray{T,2}}()

    reference_idtrend = Dict{String,CuArray{T,2}}()
    surrogate_idtrend = Dict{String,CuArray{T,2}}()

    significance = Dict{String,CuArray{T,2}}()

    for (ti_function, label) in zip(ti_function_list, ti_labels)

        reference_ti[label] = slide_estimator(X, pindctr, ti_function, window)
        surrogate_ti[label] = slide_estimator(S, pindctr, ti_function, window)
        reference_idtrend[label] = slide_idtrend(
            reference_ti[label],
            tindctr,
            pidtrend,
            trend_measure_function,
            window;
            kwargs...,
        )
        surrogate_idtrend[label] = slide_idtrend(
            surrogate_ti[label],
            tindctr,
            pidtrend,
            trend_measure_function,
            window;
            kwargs...,
        )
        significance[label] = get_percentile(
            reference_idtrend[label],
            surrogate_idtrend[label],
            ns,
            nx,
        )

    end

    result = TIresults(
        tindctr,
        tidtrend,
        S,
        reference_ti,
        surrogate_ti,
        reference_idtrend,
        surrogate_idtrend,
        significance,
    )

    return result
end

struct TIresults
    tindctr::Vector{T}
    tidtrend::Vector{T}
    S::CuArray{T,2}
    reference_ti::Dict{String,CuArray{T,2}}
    surrogate_ti::Dict{String,CuArray{T,2}}
    reference_idtrend::Dict{String,CuArray{T,2}}
    surrogate_idtrend::Dict{String,CuArray{T,2}}
    significance::Dict{String,CuArray{T,2}}
end

# TODO replace indctr by idntfr?