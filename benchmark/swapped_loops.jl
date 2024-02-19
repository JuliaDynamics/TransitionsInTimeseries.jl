using TransitionsInTimeseries

n = 10_000
t = 1:n
x = rand(n)
indicators = (var, ar1_whitenoise)
change_metrics = RidgeRegressionSlope()
config = SlidingWindowConfig(indicators, change_metrics, stride_ind = 2, stride_cha = 10)
results = estimate_changes(config, x, t)
signif = SurrogatesSignificance()
@btime significant_transitions($results, $signif)
#=
New code = swaploops-indicators-surrogates.
Old code = main branch (2023/10/30).

New code: 2.103 s (60,345 allocations: 1.12 GiB)
Old code: 2.413 s (80,415 allocations: 2.24 GiB)
=#



nseg = 10
tseg_start = 100 .+ range(0, 0.9 * n, length = nseg)
tseg_end = 900 .+ range(0, 0.9 * n, length = nseg)
config = SegmentedWindowConfig(indicators, change_metrics, collect(tseg_start),
    collect(tseg_end), min_width_cha = 50)
results = estimate_changes(config, x, t)
flags = significant_transitions(results, signif)
@btime significant_transitions($results, $signif)
#=
New code = swaploops-indicators-surrogates.
Old code = main branch (2023/10/30).

New code: 3.024 s (202,536 allocations: 956.14 MiB)
Old code: 3.881 s (403,763 allocations: 1.87 GiB)
=#