# Introduction

**TransitionIndicators.jl** is a julia package offering scalable and flexible computation of transition indicators. The latter are blueprints of critical slowing down observed before the critical transition of a multistable dynamical system.

To learn how to use this library please see [Getting started] below, and subsequently, the [Contents] page to get an overview of all offered functionality of **TransitionIndicators.jl**.

!!! tip "Star us on GitHub!"
    If you have found this library useful, please consider starring it on [GitHub](https://github.com/JuliaDynamics/TransitionIdentifiers.jl).
    This gives us an accurate lower bound of the (satisfied) user count.

## Getting started

For now, `TransitionIndicators.jl` is in its early development and can be downloaded from [here](https://github.com/JuliaDynamics/TransitionIdentifiers.jl).

## Critical slowdown

When the state of a multi-stable dynamic system experiences a change of basin of attraction, abrupt changes in the dynamics take place. For real-world systems, indicators to foresee such transitions are a valuable asset to anticipate them and mitigate their consequences. To provide such indicators, techniques of time series analysis were developed to estimate the decrease in resilience prior to a critical transition. These techniques are briefly outlined below.

In real-world applications, the forcing of a system is often subject to noise. The latter represents a disturbance, of which the decay over time delivers information on the resilience of the system for a given state. From a data analysis perspective, slower decay of disturbances lead to an increase of the:
- variance
- skewness
- kurtosis
- auto-correlation
- AR1 regression coefficient
- low-frequency power spectrum.

This is observed in the time series whenever a system is approaching a transition. This effect is known as _critical slowing down_ (CSD). The metrics listed above are here called _transition indicators_ (TIs). Other names for them can be found in literature, such as _early warning signals_, _regime-shift identifiers_ or _resilience indicators_.

## Our goals

Whereas the number of publication on CSD has drastically increased over the last years, there are relatively few open-source softwares providing high-level functions to compute TIs. The _Early Warning Signals Toolbox_ of `R` is an example of such software.

The computation of TIs is increasingly applied to fields resulting from large measurement arrays, extensive reanalysis data or verbose simulation output. This result in a multi-dimensional problem and the computational cost can therefore increase exponentially. A central aspect of the present package is to focus on performance but also on versatility. Therefore, depending on the available hardware and the problem size, the user can either perform the computation of TIs on a CPU or on a GPU while using a strictly identical syntax.

## API & terminology
