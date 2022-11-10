# [API reference](@id api_ref)

The following API reference is structured along the chronology of a typical computation.

## Utils

```@docs
structured_residual
flatten_residual
reshape_residual
get_placeholder_spacedims
```

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
    slide
    slide_estimator
    grow_window
```

### Smoothing and de-trending

```@docs
    gettrend_rollkernel
    detrend
```

### Masking

```@docs
    window_mask
    strided_window_mask
```

## [Transient indicators](@id transient_api_indicators)

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
    arp_whitenoise
    hurst_exponent
```

### Analytic regression uneven time spacing

### Numerical regression

### Frequency spectrum

```@docs
    lfps
```

## [Spatial indicators](@id transient_api_indicators)

```@docs
    spatial_variance
    spatial_skw
    spatial_krt
    eigencovar_abs
    eigencovar_rel
```


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

```@docs
    generate_stacked_fourier_surrogates
    percentile_significance
    count_positive_indicators
```

## Fast-forward way

```@docs
    indicate_transition
```