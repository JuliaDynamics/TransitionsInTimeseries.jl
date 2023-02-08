module TransitionIndicators

using LinearAlgebra
using Reexport
@reexport using StatsBase
@reexport using TimeseriesSurrogates

include("windowing.jl")
include("indicators.jl")
include("indicator_metrics.jl")

# windowing.jl
export WindowViewer

# indicators.jl
export ar1_whitenoise

# indicator_metrics.jl
export IndicatorEvolutionResults
export compute_trend
export ridge_regression
export ridge_regression_slope

end # module TransitionIndicators