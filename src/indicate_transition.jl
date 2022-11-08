"""@docs

    indicate_transition(
        ttrend::Vector{T},
        X::CuArray{T,2},
        pindctr::WindowingParams,
        pidtrend::WindowingParams,
        ns::Int,
        ti::Vector{Function},
        trend_measure::Function = ridge_regression_slope,
        window::Function = centered_wndw;
        kwargs...,
    )

Compute various transition indicators (provided in `tis`) for a residual `X` over time `ttrend`.
To track the process of transition prediction, one might want to access a verbose struct containing the surrogate data.
However, the latter can be substantial for large number of variables.
To avoid memory overflows one can choose the sparse output.
"""
function indicate_transition(
    ttrend::Vector{T},
    X::A,
    pindctr::WindowingParams,
    pidtrend::WindowingParams,
    ns::Int,
    tis::Vector{Function};
    verbose::Bool = true,
    trend_measure::Function = ridge_regression_slope,
    window::Function = centered_wndw,
    min_num_indicators::Int = length(tis),
    p_min::Real = 0.95,
    kwargs...,
) where {T<:Real, A<:Union{Matrix{T},CuArray{T,2}}}

    println(A)
    nx, nt = size(X)
    ni = length(tis)
    S = generate_stacked_fourier_surrogates(X, ns)
    ti_labels = string.(tis)

    tindctr = trim_wndw(ttrend, pindctr, window)
    tidtrend = trim_wndw(tindctr, pidtrend, window)

    reference_ti = Dict{String, A}()
    surrogate_ti = Dict{String, A}()

    reference_idtrend = Dict{String, A}()
    surrogate_idtrend = Dict{String, A}()

    significance = Dict{String, A}()

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
    positive_indicators = count_positive_indicators(indicator_trend_sigificance3D, plevel=p_min)
    positive_idx = Array(positive_indicators .>= min_num_indicators / length(tis))
    predictor_time = [tidtrend[positive_idx[i, :]] for i = 1:nx]

    if verbose
        result = verboseTIresults{T, A}(
            tindctr,
            tidtrend,
            S,
            reference_ti,
            surrogate_ti,
            reference_idtrend,
            surrogate_idtrend,
            significance,
            positive_indicators,
            predictor_time,
        )
    else
        result = sparseTIresults{T, A}(
            tindctr,
            tidtrend,
            reference_ti,
            reference_idtrend,
            significance,
            positive_indicators,
            predictor_time,
        )
    end

    return result
end
struct verboseTIresults{T<:Real, A<:Union{Matrix{T}, CuArray{T,2}}}
    tindctr::Vector{T}
    tidtrend::Vector{T}
    S::A
    reference_ti::Dict{String, A}
    surrogate_ti::Dict{String, A}
    reference_idtrend::Dict{String, A}
    surrogate_idtrend::Dict{String, A}
    significance::Dict{String, A}
    positive_indicators::A
    predictor_time::Vector{Vector{T}}
end

struct sparseTIresults{T<:Real, A<:Union{Matrix{T}, CuArray{T,2}}}
    tindctr::Vector{T}
    tidtrend::Vector{T}
    reference_ti::Dict{String, A}
    reference_idtrend::Dict{String, A}
    significance::Dict{String, A}
    positive_indicators::A
    predictor_time::Vector{Vector{T}}
end

# TODO test the new fast-forward function typing.