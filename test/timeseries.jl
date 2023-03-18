using TransitionIndicators, Test

@testset "Spacing" begin
    dt = rand()
    x_equispaced = range( 0.0, step = dt, stop = 100)
    x_nonequispaced = cumsum( rand(1000) )
    @test isequispaced(x_equispaced)
    @test !isequispaced(x_nonequispaced)
    @test equispaced_step(x_equispaced) ≈ dt
end

# Why not use this?
function my_equispaced_step(x::AbstractVector)
    return mean(diff(x))
end

x_equispaced = collect( range( 0.0, step = rand(), stop = 100) )
@btime equispaced_step(x_equispaced)
@btime my_equispaced_step(x_equispaced)