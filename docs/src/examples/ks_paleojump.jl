# # Kolmogorov-Smirnov test for detecting transitions in paleoclimate timeseries

# The goal of this example is to show how simple it is to re-create an analysis _similar_
# to what was done in the paper
# "Automatic detection of abrupt transitions in paleoclimate records",
# [Bagniewski2021](@cite). The same analysis was then used to create a database
# of transitions in paleoclimate records in [Bagniewski2023](@cite)
# Using TransitionsInTimeseries.jl and HypothesisTests.jl,
# the analysis becomes a 10-lines-of-code script (for a given timeseries).

# ## Scientific background

# The approach of [Bagniewski2021](@cite) is based on the [two-sample
# Kolmogorov Smirnov test](https://en.wikipedia.org/wiki/Kolmogorov%E2%80%93Smirnov_test#Two-sample_Kolmogorov%E2%80%93Smirnov_test).
# It tests whether the samples from two datasets or timeseries
# are distributed according to the same cumulative density function or not.
# This can be estimated by comparing the value of the _KS-statistic_
# versus some threshold that depends on the required confidence.

# The application of this test for identifying transitions in timeseries is simple:

# 1. A sliding window analysis is performed in the timesries
# 1. In each window, the KS statistic is estimated between the first half and the second
#    half of the timeseries within this window.
# 1. Transitions are defined by when the KS statistic exceeds a particular value
#    based on some confidence. The transition occurs in the middle of the window.

# We should point out that in [Bagniewski2021](@cite) the authors did a
# more detailed analysis: analyzed many different window widths, added a
# conditional clause to exclude transitions that do not exceed a predefined minimum
# "jump" in the data, and also added another conditional clause that
# filtered out transitions that are grouped in time (which is a natural consequence
# of using the Kolmogorov-Smirov test for detecting transitions).
#
# Here we won't do that post processing, mainly because it is rather simple
# to include these additional conditional clauses to filter transitions after they are found.

# ## Steps for TransitionsInTimeseries.jl

# Doing this kind of work with TransitionsInTimeseries.jl is so easy you won't even trip!
# This analysis follows the same sliding window approach showcased in our [Tutorial](@ref),
# and it even excludes the "indicator" aspect: the change metric is estimated directly
# from the input data!

# As such, we really only need to do these things before we have finished the analysis:

# 0. Obtain the input timeseries
# 1. The window size we will utilize
# 2. The input data (we will use the same example as the NGRIP data of
#    [Dansgaard-Oescher events and Critical Slowing Down](@ref) example).
# 3. Define the function that estimates the change metric (i.e., the KS-statistic)
# 3. Perform the sliding window analysis as in the [Tutorial](@ref) with [`estimate_indicator_changes`](@ref)
# 4. Estimate the "confident" transitions in the data by comparing the estimated
#    KS-statistic with a predefined threshold.

# ## Load timeseries and window length
# Following the Dansgaard-Oescher events example, we load
# the data after all the processing steps done in that example:

using DelimitedFiles

tmp = Base.download("https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/"*
    "master/timeseries/NGRIP_processed.csv")
data = readdlm(tmp)
t, xtrend, xresid, xloess = collect.(eachcol(data))

fig, ax = lines(t, xtrend; axis = (ylabel = "NGRIP (processed)", xlabel = "time"))
lines!(ax, t, xloess; linewidth = 2)
fig

# For the window, since we are using a sliding window here, we will be using a
# window of length 500 (which is approximately 1/2 to 1/4 the span between typical
# transitions found by [rasmussen-stratigraphic-2014](@cite)).

window = 500

# ## Defining the change metric function

# HypothesisTest.jl implements the Kolmogorov-Smirnov test, however here we are interested
# in the value of the test iself (the so-called KS-statistic), rather than a p-value.
# To this end, we define the following function to compute the statistic,
# which also normalizes it as in [Bagniewski2021](@cite).

using HypothesisTests

function normalized_KS_statistic(timeseries)
    N = length(timeseries)
    i = N÷2
    x = view(timeseries, 1:i)
    y = view(timeseries, (i+1):N)
    kstest = ApproximateTwoSampleKSTest(x, y)
    nx = ny = i # length of each timeseries half of total
    n = nx*ny/(nx + ny) # written fully for concreteness
    D_KS = kstest.δ # can be compared directly with sqrt(-log(α/2)/2)
    @show D_KS
    ## Rescale according to eq. (5) of the paper
    rescaled = 1 - ((1 - D_KS)/(1 - sqrt(1/n)))
    return rescaled
end

x = randn(1000)
y = 0.8randn(1000) .+ 1.0
z = randn(1000)
w = 0.6randn(1000) .- 2.0

fig, ax = density(x; color = "black")
ax.title = "showcase of normalized KS-statistic"
for q in (y, z, w)
    D_KS = normalized_KS_statistic(vcat(x, q))
    density!(ax, q; label = "D_KS = $(D_KS)")
end
axislegend(ax)
fig