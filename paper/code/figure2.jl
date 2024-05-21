using CairoMakie
using DelimitedFiles
using Random
using TransitionsInTimeseries

coefficient_of_variation(x) = std(x) / mean(x)

function main()
    # Run ewstools-tuto-1.csv first to generate ewstools_ricker.csv and ewstools_perfo.csv
    x = readdlm("../data/ewstools_ricker.csv", ',')[:, end]
    x = x[isnan.(x) .== 0]
    t = eachindex(x)

    # Choose the indicators and how to measure their change over time
    indicators = (var, coefficient_of_variation, skewness, kurtosis,
        ar1_whitenoise, LowfreqPowerSpectrum())
    stride = [1, 1, 1, 1, 1, 40]
    m = length(indicators)
    t_elapsed = zeros(m+2)

    for (i, ind) in enumerate(indicators)
        # Build configuration with adequate parameters of the sliding window
        config = SegmentedWindowConfig((ind, ind), (nothing, nothing), [t[1]], [t[end]];
            width_ind = length(x) ÷ 2, stride_ind = stride[i], whichtime = last,
            min_width_cha = 1)

        t_elapsed[i] = @belapsed estimate_changes($config, $x, $t)
    end

    config = SegmentedWindowConfig((nothing, nothing), (kendalltau, kendalltau), [t[1]], [t[end]];
        width_ind = length(x) ÷ 2, stride_ind = 1, whichtime = last, min_width_cha = 1)
    t_elapsed[m+1] = @belapsed estimate_changes($config, $x, $t)

    sgen = surrogenerator(x, BlockShuffle(), Xoshiro(1995))
    t_elapsed[m+2] = @belapsed $sgen()

    return t_elapsed
end

t_transitionsintimeries = main()
t_ewstools = readdlm("../data/ewstools_perfo.csv", ',')[:, end]
inds = eachindex(t_ewstools)
w = 0.4

fig = Figure(resolution = (600, 400))
ax = Axis(fig[1, 1])
ylims!(ax, 1e-7, 0.1)
barplot!(ax, inds .- 0.5*w, t_ewstools, label = L"ewstools $\,$", width = w,
    fillto = 1e-8)
barplot!(ax, inds .+ 0.5*w, t_transitionsintimeries, label = L"TransitionsInTimeseries.jl $\,$",
    width = w, fillto = 1e-8)
ax.yscale = log10
ax.xticks = (1:8, [L"Variance $\,$", L"Coeff. of variation $\,$", L"Skewness $\,$",
    L"Kurtosis $\,$", L"Lag-1 autocorr. $\,$", L"Spectral $\,$",
    L"Kendall $\tau$ corr. coeff.", L"Block bootstrap $\,$"])
ax.ylabel = L"Run time (s) on Ricker model data $\,$"
ax.yticks = (10.0 .^ (-8:1), [L"10^{%$e}" for e in -8:1])
ax.xgridvisible = false
ax.ygridvisible = false
ax.xticklabelrotation = π / 4
axislegend(ax, position = :lt)
save("../figures/figure2.png", fig)