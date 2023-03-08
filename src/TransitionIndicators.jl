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

include("windowing.jl")
include("indicators.jl")
include("metaanalysis.jl")
include("metaanalysis_trend.jl")
include("significance.jl")
include("load_data.jl")

# windowing.jl
export WindowViewer

# indicators.jl
export ar1_whitenoise
export std, variance, skw, krt, mean
# TODO: add lfps, restoring rate

# metaanalysis.jl
export IndicatorEvolutionResults
export SignificanceHyperParams
export analyze_indicators
export indicator_evolution
export mapwindow

# Trend metaanalysis metrics
export kendalltau, spearman     # from StatsBase
export ridge, ridge_slope
export precomputed_ridge_slope

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