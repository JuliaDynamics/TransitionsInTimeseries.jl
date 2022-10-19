module TransitionIdentifiers

using BenchmarkTools, CUDA, FFTW, StatsBase, LinearAlgebra

include("utils.jl")
include("signal_processing.jl")
include("smoothing_kernels.jl")
include("statistical_estimators.jl")
include("significance.jl")
include("highlevel_TI_computation.jl")
include("benchmark.jl")

export get_windowing_params
export centered_wndw, left_wndw, right_wndw, trim_wndw
export slide_estimator
export gettrend_rollmean, gettrend_rollkernel, detrend

export scaled_kernel

export mean, var, skw, krt
export ar1_whitenoise
export lfps

export generate_stacked_fourier_surrogates
export ridge_regression, ridge_regression_slope
export percentile_significance
export slide_idtrend

export compute_TIsignificance

export benchmark_functions_over_size
export benchmark_cpu_vs_gpu
export lineplot_benchmark

end # module TransitionIdentifiers
