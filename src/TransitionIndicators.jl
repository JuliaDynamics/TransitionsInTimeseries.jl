module TransitionIndicators

using LinearAlgebra
using Reexport
using Random
using StatsBase
@reexport using TimeseriesSurrogates

include("windowing.jl")
include("indicators.jl")
include("evolution_metrics.jl")
include("metaanalysis_trend.jl")
include("significance.jl")

# windowing.jl
export WindowViewer

# indicators.jl
export ar1_whitenoise

# evolution_metrics.jl
export IndicatorEvolutionResults
export MetaAnalysisParameters
export init_metaanalysis_params
export analyze_indicators
export analyze_indicator

# Trend metaanalysis metrics
export corspearman, corkendall # from StatsBase
export ridge, ridge_slope

# significance.jl
export measure_significance
export measure_significances
export normalized_confidence_intervall
export which_percentile
export percentile_idx
export normalized_percentile_distance
export intround

end # module TransitionIndicators