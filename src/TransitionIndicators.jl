module TransitionIndicators

using LinearAlgebra
using Reexport
@reexport using StatsBase
@reexport using TimeseriesSurrogates

include("windowing.jl")
include("indicators.jl")
include("evolution_metrics.jl")
include("significance.jl")

# windowing.jl
export WindowViewer

# indicators.jl
export ar1_whitenoise

# evolution_metrics.jl
export IndicatorEvolutionResults
export indicator_evolution

export ridge
export ridge_slope
export precompute_ridge
export precompute_ridge_slope
export precomputed_ridge_slope

# significance.jl
export measure_significance
export gaussian_quantile
export which_quantile
export quantile_idx
export normalized_quantile_distance
export intround

end # module TransitionIndicators