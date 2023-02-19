# TransitionIndicators.jl

```@docs
TransitionIndicators
```

## Example
```@example MAIN
using TransitionIndicators
using Statistics: mean, var

n = 1001
t = collect(1.0:n)
x = copy(t)
p = init_metaanalysis_params(n_surrogates = 100)
res = analyze_indicators(t, x, [mean, var], ridge_slope, p)

# The trend of mean(windowview) is the stride for x=t
meantrend_ground_truth = fill(p.wv_indicator_stride, length(res.t_evolution))
# The trend of var(windowview) is 0 for x any affine function of t.
vartrend_ground_truth = fill(0.0, length(res.t_evolution))
```


## API

```@docs
analyze_indicators
```
