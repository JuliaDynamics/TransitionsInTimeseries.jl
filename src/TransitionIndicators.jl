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
using DelimitedFiles

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
# TODO: permutation entropy
export kendalltau, spearman     # from StatsBase
export RidgeRegression

# analysis
export IndicatorsConfig, SignificanceConfig, IndicatorsResults
export indicators_analysis, indicators_significance

# timeseries
export isequispaced, equispaced_step

# significance
export significant, Quantile, ThresholdQuantile, Sigma, ThresholdSigma

# load_data.jl
export load_linear_vs_doublewell

end # module TransitionIndicators