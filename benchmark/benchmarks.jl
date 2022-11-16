using TransitionIndicators
using Statistics, BenchmarkTools
function inplace_map!(f::Function, y::Vector{T}, wv::WindowViewer) where {T<:Real}
    @inbounds for (i, windowview) in enumerate(wv)
        y[i] = f(windowview)
    end
end

n = 100_000
x = collect(1:n)
halfwidth = 2
stride = 2
wv = WindowViewer(x, halfwidth, stride)

windowed_views = collect(wv)
var_x = zeros(length(wv))

@btime map(var, $wv)
@btime map!(var, $var_x, $windowed_views)
@btime inplace_map!(var, $var_x, $wv)
