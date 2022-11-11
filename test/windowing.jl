using TransitionIndicators, Test, StatsBase

@testset "sliding variance over vector" begin
    n = 100
    x = collect(1:n)
    halfwidth = Int(ceil(rand() * (n/2-1)))     # Draw integer on [1, n/2-1]. Any of those should work.
    stride = 2
    ground_truth = var(1:2*halfwidth+1)         # x is a linear function (x = t) --> variance independent of window.

    brackets = [left, center, right]
    window_viewers = [WindowViewer(x, halfwidth, stride, bracket) for bracket in brackets]
    variance_window_viewers = vcat([map(StatsBase.var, wv)' for wv in window_viewers]...)

    # Check if all entries of result are equal to ground truth.
    @test sum( variance_window_viewers .== ground_truth ) == length(variance_window_viewers)
end
