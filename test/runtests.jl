using Test

defaultname(file) = uppercasefirst(replace(splitext(basename(file))[1], '_' => ' '))
testfile(file, testname=defaultname(file)) = @testset "$testname" begin; include(file); end

@testset "TranstionIndicators.jl" begin
    testfile("dummy.jl")
    testfile("windowing.jl")
end
