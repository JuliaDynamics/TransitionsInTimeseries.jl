module TransitionIndicators

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end TransitionIndicators

using LinearAlgebra
using Random
using StatsBase
using Downloads
using CSV
using DataFrames

using Reexport
@reexport using TimeseriesSurrogates

include("misc/windowing.jl")
include("misc/timeseries.jl")
include("misc/load_data.jl")

include("library/indicators.jl")
include("library/change_metrics_trend.jl")

include("analysis/analysis_types.jl")
include("analysis/perform_analysis.jl")
include("analysis/significance.jl")

# windowing.jl
export WindowViewer, windowmap, windowmap!, midpoint, midvalue

# library
export ar1_whitenoise
export mean, std, var, skewness, kurtosis
# TODO: add lfps, restoring rate, permutation entropy
export kendalltau, spearman     # from StatsBase
export RidgeRegression

# analysis
export IndicatorsConfig, SignificanceConfig, IndicatorsResults
export indicators_analysis, indicators_significance


# export threshold_indicators
# export measure_significance
# export confidence_interval
# export normalized_percentile
# export which_percentile
# export percentile_idx
# export intround

# load_data.jl
export load_linear_vs_doublewell

end # module TransitionIndicators