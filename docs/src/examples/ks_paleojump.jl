# # Kolmogorov-Smirnov test for detecting transitions

# The goal of this example is to show how simple it is to re-create an analysis _similar_
# to what was done in the paper
# "Automatic detection of abrupt transitions in paleoclimate records"
# by Bagniewski, Ghil, Rousseau, with DOI: <https://doi.org/10.1063/5.0062543>.
# In fact, using TransitionsInTimeseries.jl and HypothesisTests.jl,
# the whole analysis becomes a 10-lines-of-code script
# (for each timeseries in the Paleojump record).


# ## Scientific background

# The approach of the paper is based on the [two-sample
# Kolmogorov Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test#Two-sample_Kolmogorov%E2%80%93Smirnov_test).
# It tests whether the samples from two datasets
# are distributed according to the same PDF or not.
# This can be estimated by comparing the value of the _KS-statistic_
# versus some threshold that depends on the required confidence.

# The application of this test for identifying transitions in timeseries is simple:

# 1. A sliding window analysis is performed in the timesries
# 1. In each window, the KS statistic is estimated between the first half and the second
#    half of the timeseries within this window.
# 1. Transitions are defined by when the KS statistic exceeds a particular value
#    based on some confidence. The transition occurs in the middle of the window.

# We should point out that the publication we are referring to did a much
# more detailed analysis: it analyzed many different window widths, and added a
# conditional clause to exclude transitions that do not exceed a predefined minimum
# "jump" in the data. Here we won't do that (mainly because it is rather simple
# to include these additional conditional clauses to filter transitions after they are found).



# ## Defining the change metric function


# HypothesisTest.jl implements this test, however here we are interested
# in the value of the test iself (the so-called KS-statistic), rather than a p-value.
# So we will need to


using HypothesisTests

x = randn(1000)
y = 0.9randn(1000) .+ 0.1

function ks_statistic(x, y)
    kstest = ApproximateTwoSampleKSTest(x, y)
    nx, ny = length(x), length(y)
    n = nx*ny/(nx + ny)
    return sqrt(n)*kstest.δ
end