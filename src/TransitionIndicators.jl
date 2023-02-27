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

using Reexport
@reexport using TimeseriesSurrogates

include("windowing.jl")
include("indicators.jl")
include("metaanalysis.jl")
include("metaanalysis_trend.jl")
include("significance.jl")

# windowing.jl
export WindowViewer

# indicators.jl
export ar1_whitenoise

# metaanalysis.jl
export IndicatorEvolutionResults
export MetaAnalysisParameters
export init_metaanalysis_params
export analyze_indicators
export analyze_indicator
export mapwindow

# Trend metaanalysis metrics
export spearman, kendalltau # from StatsBase
export ridge, ridge_slope
export precompute_ridge_slope, precomputed_ridge_slope
export precompute_ridge

# significance.jl
export measure_significance
export measure_significances
export normalized_confidence_intervall
export which_percentile
export percentile_idx
export normalized_percentile_distance
export intround

end # module TransitionIndicators