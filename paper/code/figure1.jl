using TransitionsInTimeseries, DelimitedFiles, CairoMakie, Random

# Run ewstools-tuto-1.csv first to generate ewstools_ricker.csv
x = readdlm("../data/ewstools_ricker.csv", ',')[:, end]
x = x[isnan.(x) .== 0]
t = eachindex(x)

# Choose the indicators and how to measure their change over time
indicators = (var, ar1_whitenoise)
change_metrics = (kendalltau, kendalltau)
config = SegmentedWindowConfig(indicators, change_metrics, [t[1]], [t[end]];
    width_ind = length(x) ÷ 2, whichtime = last, min_width_cha = 50)
results = estimate_changes(config, x, t)
signif = SurrogatesSignificance(n = 1000, tail = [:right, :right], rng = Xoshiro(1995))
flags = significant_transitions(results, signif)
fig = plot_changes_significance(results, signif)
save("../figures/figure1.png", fig)