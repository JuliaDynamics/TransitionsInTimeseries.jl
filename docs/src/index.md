# TransitionIndicators.jl

<!-- ```@docs
TransitionIndicators
``` -->

**TransitionIndicators.jl** is a julia package offering scalable and flexible computation of transition indicators. The latter are blueprints of critical slowing down observed before the critical transition of a multistable dynamical system.

To learn how to use this library please see [Getting started](@ref) below, and subsequently, the [Contents](@ref) page to get an overview of all offered functionality of **TransitionIndicators.jl**.

!!! tip "Star us on GitHub!"
    If you have found this library useful, please consider starring it on [GitHub](https://github.com/JuliaDynamics/TransitionIdentifiers.jl).
    This gives us an accurate lower bound of the (satisfied) user count.

## Getting started

For now, `TransitionIndicators.jl` is in its early development and can be downloaded from (https://github.com/JuliaDynamics/TransitionIdentifiers.jl)[here].

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



### Probabilities

Entropies and other complexity measures are typically computed based on _probability distributions_.
These are obtained from [Input data for Entropies.jl](@ref) in a plethora of different ways.
The central API function that returns a probability distribution (in fact, just a vector of probabilities) is [`probabilities`](@ref), which takes in a subtype of [`ProbabilitiesEstimator`](@ref) to specify how the probabilities are computed.
All estimators available in Entropies.jl can be found in the [estimators page](@ref probabilities_estimators).

### Entropies

Entropy is an established concept in statistics, information theory, and nonlinear dynamics.
However it is also an umbrella term that may mean several computationally different quantities.
In Entropies.jl, we provide the generic function [`entropy`](@ref) that tries to both clarify the disparate "entropy concepts", while unifying them under a common interface that highlights the modular nature of the word "entropy".

Most of the time, computing an entropy boils down to two simple steps: first estimating a probability distribution, and then applying one of the so-called "generalized entropy" formulas to the distributions.
Thus, any of the implemented [probabilities estimators](@ref probabilities_estimators) can be used to compute generalized entropies.

!!! tip "There aren't many entropies, really."
    A crucial thing to clarify is that many quantities that are named as entropies (e.g., permutation entropy [`entropy_permutation`](@ref), wavelet entropy [`entropy_wavelet`](@ref), etc.), are _not really new entropies_. They are new probabilities estimators. They simply devise a new way to calculate probabilities from data, and then plug those probabilities into formal entropy formulas such as the Shannon entropy. The probabilities estimators are smartly created so that they elegantly highlight important aspects of the data relevant to complexity.

    These names are common place, and so in Entropies.jl we provide convenience functions like [`entropy_wavelet`](@ref). However, it should be noted that these functions really aren't anything more than 2-lines-of-code wrappers that call [`entropy`](@ref) with the appropriate [`ProbabilitiesEstimator`](@ref).

    There are only a few exceptions to this rule, which are quantities that are able to compute Shannon entropies via alternate means, without explicitly computing some probability distributions. These are `IndirectEntropy` instances, such as [`Kraskov`](@ref).

### Complexity measures

Other complexity measures, which strictly speaking don't compute entropies, and may or may not explicitly compute probability distributions, appear in the [Complexity measures](@ref complexity_measures) section.

## Input data for Entropies.jl

The input data type typically depend on the probability estimator chosen. In general though, the standard DynamicalSystems.jl approach is taken and as such we have three types of input data:

- _Timeseries_, which are `AbstractVector{<:Real}`, used in e.g. with [`WaveletOverlap`](@ref).
- _Multi-dimensional timeseries, or datasets, or state space sets_, which are [`Dataset`](@ref), used e.g. with [`NaiveKernel`](@ref).
- _Spatial data_, which are higher dimensional standard `Array`s, used e.g. with  [`SpatialSymbolicPermutation`](@ref).
