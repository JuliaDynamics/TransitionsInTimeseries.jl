using HypothesisTests

function kolmogorov_smirnov(x)
    m = length(x)÷2
    x1 = view(x, 1:m)
    x2 = view(x, (m+1):length(x))
    test = ApproximateTwoSampleKSTest(x1, x2)
    return pvalue(test)
end