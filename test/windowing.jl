using TransitionIndicators, Test, Statistics

@testset "sliding variance over vector" begin
    n = 100
    x = collect(1:n)
    indicator_window = (width = Int(ceil(rand() * (n/2-1))), stride = 2)
    windowed_time = windowmap(maximum, x; indicator_window...)
    windowed_variance = windowmap(Statistics.var, x; indicator_window...)

    # x is a linear function (x = t) -→ variance independent of window.
    true_time = x[indicator_window.width+1:indicator_window.stride:end]
    true_variance = var(1:indicator_window.width+1)
    @test sum( windowed_variance .≈ true_variance ) == length(windowed_variance)
    @test windowed_time == true_time
end