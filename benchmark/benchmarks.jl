using TransitionsInTimeseries
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


#################################################
# surrogates_loop!
#################################################

using TransitionsInTimeseries

n = 1001
t = collect(1.0:n)
x = copy(t)

indicators = (mean, var)
change_metrics = RidgeRegressionSlope()

config = SlidingWindowConfig(indicators, change_metrics;
    width_ind = 100, stride_ind = 1,
    width_cha = 100, stride_cha = 1, whichtime = last,
)

res = estimate_indicator_changes(config, x, t)
signif = SurrogatesSignificance(n = 100, tail = :both, p = 0.1)
flags = significant_transitions(res, signif)

@btime estimate_indicator_changes($config, $x, $t)
# 141.557 μs (4 allocations: 40.56 KiB)
@btime significant_transitions($res, $signif)
#=
Former version of the code: 5.178 ms (1954 allocations: 3.03 MiB)
Current version of the code: 4.863 ms (738 allocations: 2.96 MiB)

Although there is an improvement in the number of allocations, the performance
is not significantly improved.
=#