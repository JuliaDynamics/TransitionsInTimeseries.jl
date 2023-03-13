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

include("analysis/metaanalysis.jl")
include("analysis/significance.jl")

# windowing.jl
export WindowViewer, windowmap

# indicators.jl
export ar1_whitenoise, midpoint
export mean, std, var, skewness, kurtosis
# TODO: add lfps, restoring rate, permutation entropy

# metaanalysis.jl
export IndicatorEvolutionResults
export SignificanceHyperParams
export analyze_indicators
export indicator_evolution

# Trend metaanalysis metrics
export kendalltau, spearman     # from StatsBase
export RidgeRegression

# significance.jl
export threshold_indicators
export measure_significance
export confidence_interval
export normalized_percentile
export which_percentile
export percentile_idx
export intround

# load_data.jl
export load_linear_vs_doublewell

end # module TransitionIndicators