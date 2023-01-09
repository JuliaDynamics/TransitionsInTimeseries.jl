using TransitionIndicators
using Statistics, BenchmarkTools

#################################################
# Windowing
#################################################

function inplace_map!(f::Function, y::Vector{T}, wv::WindowViewer) where {T<:Real}
    @inbounds for (i, windowview) in enumerate(wv)
        y[i] = f(windowview)
    end
end

n = 100_000
x = collect(1:n)
width = 5
stride = 2
wv = WindowViewer(x, width, stride)

windowed_views = collect(wv)
var_x = zeros(length(wv))

@btime map(var, $wv)
@btime map!(var, $var_x, $windowed_views)
@btime inplace_map!(var, $var_x, $wv)
#=
Output on Jan's machine:
  1.046 ms (2 allocations: 390.67 KiB)
  969.439 μs (0 allocations: 0 bytes)
  1.148 ms (0 allocations: 0 bytes)
=#

#################################################
# AR1
#################################################

# Check 0 allocations
@btime inplace_map!(ar1_whitenoise, $var_x, $wv)
