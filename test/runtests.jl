using Test

"""

    load_linear_vs_doublewell()

Load prototypical data from a linear and a double-well model to test some indicators.
"""
function load_linear_vs_doublewell()
    url = "https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/"*
        "timeseries/linear_vs_doublewell.csv"
    tmp = Downloads.download(url)
    data = readdlm(tmp, ',', Float64, skipstart = 1)
    return view(data, :, 1), view(data, :, 2), view(data, :, 3) # t, xlin, xnlin
end

defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))
testfile(file, testname=defaultname(file)) = @testset "$testname" begin; include(file); end

@testset "TranstionIndicators.jl" begin
    testfile("timeseries.jl")
    testfile("windowing.jl")
    testfile("indicators.jl")
    testfile("change_metrics.jl")
    testfile("full_analysis.jl")
    testfile("slope_change.jl")
end
