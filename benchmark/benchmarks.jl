using StatsBase, TransitionIndicators

# Needs to be computed only once for a given timeseries.
# Embed in src?
function vectorize_windowviewer(wv::WindowViewer)
    return [windowview for windowview in wv]
end

begin
    n = 1_000_000
    x = collect(1:n)
    halfwidth = 2                               # Draw integer on [1, n/2-1]. Any of those should work.
    stride = 2

    wv = WindowViewer(x, halfwidth, stride)
    vectorized_window = vectorize_windowviewer(wv)
    var_x = zeros(length(wv.strided_indices))

    @btime map($StatsBase.var, $wv)
    @btime map!($StatsBase.var, $var_x, $vectorized_window)
end

