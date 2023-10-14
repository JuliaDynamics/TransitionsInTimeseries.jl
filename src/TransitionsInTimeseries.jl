module TransitionsInTimeseries

# Use the README as the module docs
@doc let
    path = joinpath(dirname(@__DIR__), "README.md")
    include_dependency(path)
    read(path, String)
end TransitionsInTimeseries

using LinearAlgebra
using Random
using Downloads
using DelimitedFiles
using InteractiveUtils
using FFTW

using Reexport
@reexport using TimeseriesSurrogates

include("misc/params.jl")
include("misc/windowing.jl")
include("misc/timeseries.jl")
include("misc/load_data.jl")
include("misc/precomputation.jl")

include("analysis/api.jl")
include("analysis/sliding_window.jl")
include("analysis/segmented_window.jl")
include("significance/api_significance.jl")
include("significance/surrogates_significance.jl")
include("significance/quantile_significance.jl")

include("indicators/critical_slowing_down.jl")
include("indicators/distribution_distance.jl")
include("indicators/nlts.jl")
include("indicators/spectral.jl")
include("indicators/statistics.jl")

include("change_metrics/slope.jl")
include("change_metrics/valuediff.jl")


# windowing.jl
export WindowViewer, windowmap, windowmap!, midpoint, midvalue

# library
export PrecomputableFunction, precompute
export ar1_whitenoise
export LowfreqPowerSpectrum, PrecomputedLowfreqPowerSpectrum
export mean, std, var, skewness, kurtosis # from StatsBase
export permutation_entropy
export kendalltau, spearman
export ridgematrix, RidgeRegressionSlope, PrecomputedRidgeRegressionSlope
export difference_of_means

# analysis
export IndicatorsChangesConfig, SlidingWindowConfig, SegmentedWindowConfig
export SlidingWindowResults, SegmentedWindowResults
export estimate_indicator_changes, IndicatorsChangesResults
export TransitionsSignificance, significant_transitions, segmented_significance
export QuantileSignificance, SigmaSignificance, SurrogatesSignificance

# timeseries
export isequispaced, equispaced_step

# load_data.jl
export load_linear_vs_doublewell

end # module TransitionsInTimeseries