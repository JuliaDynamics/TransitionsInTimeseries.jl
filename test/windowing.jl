using TransitionIndicators, Test, StatsBase

@testset "sliding variance over vector" begin
    n = 100
    x = collect(1:n)
    halfwidth = Int(ceil(rand() * (n/2-1)))     # Draw integer on [1, n/2-1]. Any of those should work.
    stride = 2
    ground_truth = var(1:2*halfwidth+1)         # x is a linear function (x = t) --> variance independent of window.

    wv = WindowViewer(x, halfwidth, stride)
    variance_window_viewer = map(StatsBase.var, wv)

    # Check if all entries of result are equal to ground truth.
    @test sum( variance_window_viewer .== ground_truth ) == length(variance_window_viewer)
end
