module TransitionIdentifiers

# TODO: reduce to only needed functions
using BenchmarkTools
using CUDA
using FFTW
using StatsBase
using LinearAlgebra

include("utils.jl")
include("signal_processing.jl")
include("smoothing_kernels.jl")
include("indicators.jl")
include("trend_estimation.jl")
include("significance.jl")
include("indicate_transition.jl")
include("benchmark.jl")

export get_windowing_params
export centered_wndw, left_wndw, right_wndw, trim_wndw
export slide_estimator
export gettrend_rollmean, get_trend, detrend
export strided_window_mask

export uniform_kernel, triangular_kernel, parabolic_kernel
export biweight_kernel, triweight_kernel, tricube_kernel
export gaussian_kernel, cosine_kernel, logistic_kernel
export sigmoid_kernel, scaled_kernel

export mean, masked_mean, var, masked_meansquare, skw, krt
export ar1_whitenoise, masked_ar1_whitenoise
export lfps

export generate_stacked_fourier_surrogates
export ridge_regression, ridge_regression_slope
export kendall_tau, scaled_kendall_tau
export percentile_significance
export slide_idtrend
export count_positive_indicators

export stack_indicators, indicate_transition

export benchmark_functions_over_size
export benchmark_cpu_vs_gpu
export lineplot_benchmark

end # module TransitionIdentifiers
