"""
    function TIs_significance(
        t::Vector{T},
        X::CuArray{T,2},
        pindctr::WindowingParams,
        pidtrend::WindowingParams,
        ns::Int,
        ti::Vector{Function},
        trend_measure::Function = ridge_regression_slope,
        window::Function = centered_wndw;
        kwargs...,
    )

Compute various transition identifiers (provided in `tis`) for a residual `X` over time `t`.
To track the process of transition prediction, one might want to access a verbose struct containing the surrogate data.
However, the latter can be substantial for large number of variables.
To avoid memory overflows one can choose the sparse output.
"""

function identify_transition(
    t::Vector{T},
    X::Union{Matrix{T}, CuArray{T,2}},
    pindctr::WindowingParams,
    pidtrend::WindowingParams,
    ns::Int,
    tis::Vector{Function};
    verbose::Bool = true,
    trend_measure::Function = ridge_regression_slope,
    window::Function = centered_wndw,
    kwargs...,
) where {T}

    nx, nt = size(X)
    ni = length(tis)
    S = generate_stacked_fourier_surrogates(X, ns)
    ti_labels = string.(tis)

    tindctr = trim_wndw(t, pindctr, window)
    tidtrend = trim_wndw(tindctr, pidtrend, window)

    reference_ti = Dict{String,CuArray{T,2}}()
    surrogate_ti = Dict{String,CuArray{T,2}}()

    reference_idtrend = Dict{String,CuArray{T,2}}()
    surrogate_idtrend = Dict{String,CuArray{T,2}}()

    significance = Dict{String,CuArray{T,2}}()

    for (ti_function, label) in zip(tis, ti_labels)

        reference_ti[label] = slide_estimator(X, pindctr, ti_function, window)
        surrogate_ti[label] = slide_estimator(S, pindctr, ti_function, window)
        reference_idtrend[label] = slide_idtrend(
            reference_ti[label],
            tindctr,
            pidtrend,
            trend_measure,
            window;
            kwargs...,
        )
        surrogate_idtrend[label] = slide_idtrend(
            surrogate_ti[label],
            tindctr,
            pidtrend,
            trend_measure,
            window;
            kwargs...,
        )
        significance[label] = percentile_significance(
            reference_idtrend[label],
            surrogate_idtrend[label],
            ns,
            nx,
        )

    end

    indicator_trend_sigificance = [significance[label] for label in ti_labels]
    indicator_trend_sigificance3D = stack_indicators(indicator_trend_sigificance)
    predictor_idx = predict_transition(indicator_trend_sigificance3D, nindicators=ni)
    predictor_time = [tidtrend[predictor_idx[i,:]] for i in 1:nx]

    if verbose
        result = verboseTIresults{T}(
            tindctr,
            tidtrend,
            S,
            reference_ti,
            surrogate_ti,
            reference_idtrend,
            surrogate_idtrend,
            significance,
            predictor_time,
        )
    else
        result = sparseTIresults{T}(
            tindctr,
            tidtrend,
            reference_ti,
            reference_idtrend,
            significance,
            predictor_time,
        )
    end

    return result
end
struct verboseTIresults{T}
    tindctr::Vector{T}
    tidtrend::Vector{T}
    S::Union{Matrix{T}, CuArray{T,2}}
    reference_ti::Dict{String, Union{Matrix{T}, CuArray{T,2}}}
    surrogate_ti::Dict{String, Union{Matrix{T}, CuArray{T,2}}}
    reference_idtrend::Dict{String, Union{Matrix{T}, CuArray{T,2}}}
    surrogate_idtrend::Dict{String, Union{Matrix{T}, CuArray{T,2}}}
    significance::Dict{String, Union{Matrix{T}, CuArray{T,2}}}
    predictor::Vector{Vector{T}}
end

struct sparseTIresults{T}
    tindctr::Vector{T}
    tidtrend::Vector{T}
    reference_ti::Dict{String, Union{Matrix{T}, CuArray{T,2}}}
    reference_idtrend::Dict{String, Union{Matrix{T}, CuArray{T,2}}}
    significance::Dict{String, Union{Matrix{T}, CuArray{T,2}}}
    predictor::Vector{Vector{T}}
end

# TODO replace indctr by idntfr?