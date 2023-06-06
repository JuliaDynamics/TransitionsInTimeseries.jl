using TransitionsInTimeseries, Test, Statistics

@testset "sliding variance over vector" begin
    n = 1_000
    x = collect(1:n)
    w, s = 23, 2
    windowed_time = windowmap(last, x; width = w, stride = s)
    windowed_variance = windowmap(Statistics.var, x; width = w, stride = s)

    # x is a linear function (x = t) -→ variance independent of window.
    true_time = x[w:s:end]
    true_variance = fill(var(1:w), length(windowed_time))
    @test windowed_variance == true_variance
    @test windowed_time == true_time
end