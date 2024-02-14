using TransitionsInTimeseries, DelimitedFiles, CairoMakie, Random

coefficient_of_variation(x) = std(x) / mean(x)

function main()
    x = readdlm("ewstools-tuto-1.csv", ',')[:, end]
    t = eachindex(x)

    # Choose the indicators and how to measure their change over time
    indicators = (var, coefficient_of_variation, skewness, kurtosis,
        ar1_whitenoise, LowfreqPowerSpectrum())
    stride = [1, 1, 1, 1, 1, 40]
    n, m = 100, length(indicators)
    t_elapsed = zeros(m+2)


    for (i, ind) in enumerate(indicators)
        # Build configuration with adequate parameters of the sliding window
        config = SegmentedWindowConfig((ind, ind), (nothing, nothing), [t[1]], [t[end]];
            width_ind = length(x) ÷ 2, stride_ind = stride[i], whichtime = last,
            min_width_cha = 1)

        t0 = time()
        for i in 1:n
            # Compute the metrics over sliding windows and their significance
            results = estimate_indicator_changes(config, x, t)
        end
        t_elapsed[i] = (time() - t0) / 2
    end

    config = SegmentedWindowConfig((nothing, nothing), (kendalltau, kendalltau), [t[1]], [t[end]];
        width_ind = length(x) ÷ 2, stride_ind = 1, whichtime = last, min_width_cha = 1)
    t0 = time()
    for i in 1:n
        results = estimate_indicator_changes(config, x, t)
    end
    t_elapsed[m+1] = (time() - t0) / 2

    sgen = surrogenerator(x, BlockShuffle(), Xoshiro(1995))
    t0 = time()
    for i in 1:n
        s = sgen()
    end
    t_elapsed[m+2] = time() - t0

    return t_elapsed
end

t_tt = main()
t_et = [0.03840542, 0.05554581, 0.03895116, 0.04029274, 7.96556187,
    2.73067856, 0.39529872, 0.02751493]

# [0.04681492, 8.13679838, 0.04035759, 0.09219241]
inds = eachindex(t_et)
w = 0.4

fig, ax = barplot(inds .- 0.5*w, t_et, label = L"ewstools $\,$", width = w,
    fillto = 1e-5)
barplot!(ax, inds .+ 0.5*w, t_tt, label = L"TransitionsInTimeseries.jl $\,$",
    width = w, fillto = 1e-5)
ax.yscale = log10
ax.xticks = (1:8, [L"Variance $\,$", L"Coeff. of variation $\,$", L"Skewness $\,$",
    L"Kurtosis $\,$", L"Lag-1 autocorr. $\,$", L"Spectral $\,$",
    L"Kendall $\tau$ corr. coeff.", L"Block bootstrap $\,$"])
ax.ylabel = L"Run time of 100 computations (s) $\,$"
ax.yticks = (10.0 .^ (-5:1), [L"10^{%$e}" for e in -5:1])
ax.xgridvisible = false
ax.ygridvisible = false
ax.xticklabelrotation = π / 4
axislegend(ax, position = :lt)
save("../figures/figure2.png", fig)