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
using InteractiveUtils
using FFTW
using ComplexityMeasures

using Reexport
@reexport using TimeseriesSurrogates

include("misc/params.jl")
include("misc/windowing.jl")
include("misc/timeseries.jl")
include("misc/load_data.jl")

include("library/indicators.jl")
include("library/change_metrics_trend.jl")

include("analysis/analysis_types.jl")
include("analysis/perform_analysis.jl")

# params.jl
export IndicatorsParams, ChangeMetricsParams, init_metrics

# windowing.jl
export WindowViewer, windowmap, windowmap!, midpoint, midvalue

# library
export IndicatorsParams, ChangeMetricsParams
export ar1_whitenoise, LowfreqPowerSpectrum
export mean, std, var, skewness, kurtosis       # from StatsBase
export SymbolicPermutation, entropy_normalized  # from ComplexityMeasures
export PermutationEntropy
export kendalltau, spearman
export RidgeRegressionSlope

# analysis
export IndicatorsConfig, SignificanceConfig, IndicatorsResults
export indicators_analysis, indicators_significance

# timeseries
export isequispaced, equispaced_step

# load_data.jl
export load_linear_vs_doublewell

end # module TransitionIndicators