using TransitionIndicators, Test, Statistics

@testset "sliding variance over vector" begin
    n = 100
    x = collect(1:n)
    width = Int(ceil(rand() * (n-1)))     # Draw integer on [1, n/2-1]. Any of those should work.
    stride = 2
    ground_truth = var(1:width+1)         # x is a linear function (x = t) -➡ variance independent of window.

    wv = WindowViewer(x, width, stride)
    variance_window_viewer = map(Statistics.var, wv)

    # Check if all entries of result are equal to ground truth.
    @test sum( variance_window_viewer .== ground_truth ) == length(variance_window_viewer)
end
