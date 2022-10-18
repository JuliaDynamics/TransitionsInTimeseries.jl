module TransitionIdentifiers

using BenchmarkTools, CUDA, FFTW, StatsBase, LinearAlgebra

include("utils.jl")
include("signal_processing.jl")
include("smoothing_kernels.jl")
include("statistical_estimators.jl")
include("significance.jl")
include("highlevel_TI_computation.jl")
include("benchmark.jl")

export skw

end # module TransitionIdentifiers
