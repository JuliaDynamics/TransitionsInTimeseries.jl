module TransitionIdentifiers

using LinearAlgebra
using CUDA
using FFTW: rfft, irfft
using StatsBase: corkendall
using SparseArrays: spzeros, spdiagm, SparseMatrixCSC

include("utils.jl")
include("signal_processing.jl")
include("smoothing_kernels.jl")
include("indicators.jl")
include("trend_estimation.jl")
include("significance.jl")
include("indicate_transition.jl")

#############################
# Utils
#############################

export Residual
export structured_residual
export flatten_residual
export reshape_residual
export get_placeholder_spacedims

#############################
# Signal processing
#############################

export get_windowing_params
export centered_wndw, left_wndw, right_wndw
export trim_wndw
export slide_estimator
export grow_window
export gettrend_rollmean
export gettrend_rollkernel
export get_trend
export detrend
export window_mask
export strided_window_mask

#############################
# Kernels
#############################

export uniform_kernel
export triangular_kernel
export parabolic_kernel
export biweight_kernel
export triweight_kernel
export tricube_kernel
export gaussian_kernel
export cosine_kernel
export logistic_kernel
export sigmoid_kernel
export scaled_kernel

#############################
# Indicators
#############################

export mean_lastdim
export masked_mean_lastdim
export var
export masked_meansquare
export std
export skw
export krt
export ar1_whitenoise
export masked_ar1_whitenoise
export arp_whitenoise
export hurst_exponent
export lfps

#############################
# Significance
#############################

export generate_stacked_fourier_surrogates
export ridge_regression
export ridge_regression_slope
export kendall_tau
export scaled_kendall_tau
export percentile_significance
export slide
export count_positive_indicators
export stack_indicators

#############################
# High-level interface
#############################

export indicate_transition

end # module TransitionIdentifiers
