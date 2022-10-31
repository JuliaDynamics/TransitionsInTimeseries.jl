# [API reference](@id api_ref)

The following API reference is structured along the chronology of a typical computation.

## Signal processing

### Windowing

```@docs
    get_windowing_params
```

```@docs
    left_wndw
```

```@docs
    centered_wndw
```

```@docs
    right_wndw
```

```@docs
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
```

```@docs
    detrend
```




## [Indicators](@id api_indicators)

### Statistical moments

```@docs
    mean
```

```@docs
    masked_mean
```

```@docs
    var
```

```@docs
    masked_meansquare
```

```@docs
    skw
```

```@docs
    krt
```

### Analytic regression

```@docs
    ar1_whitenoise
```

```@docs
    masked_ar1_whitenoise
```

### Analytic regression uneven time spacing

### Numerical regression

### Frequency spectrum

```@docs
    lfps
```

### Spatial indicators




## Trend estimation

### Linear regression

```@docs
    ridge_regression
```

```@docs
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