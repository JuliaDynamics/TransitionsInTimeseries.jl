# [API reference](@id api_ref)

The following API reference is structured along the chronology of a typical computation.

## Signal processing

### Windowing

```@docs
    get_windowing_params
    left_wndw
    centered_wndw
    right_wndw
    trim_wndw
```

## Kernels

```@docs
    scaled_kernel
```

### Sliding function over arrays

```@docs
    slide_estimator
```

### Smoothing and de-trending

```@docs
    gettrend_rollkernel
    detrend
```




## [Indicators](@id api_indicators)

### Statistical moments

```@docs
    mean_lastdim
    masked_mean_lastdim
    var
    masked_meansquare
    skw
    krt
```

### Analytic regression


```math
\theta = \frac{\sum_i x_i \, x_{i+1}}{\sum_i x_i \, x_i}
```


```@docs
    ar1_whitenoise
    masked_ar1_whitenoise
```

### Analytic regression uneven time spacing

### Numerical regression

### Frequency spectrum

```@docs
    lfps
```

### Spatial indicators




## [Trend estimation](@id api_trends)

### Linear regression

```@docs
    ridge_regression
    ridge_regression_slope
```

### Correlation

```@docs
    kendall_tau
```





## Significance computation

### Percentile significance

```@docs
    percentile_significance
```

### Counting positive indicators

```@docs
    count_positive_indicators
```




## Fast-forward way

```@docs
    indicate_transition
```