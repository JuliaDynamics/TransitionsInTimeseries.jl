# TransitionIndicators.jl

![TransitionIndicators.jl](https://github.com/JuliaDynamics/JuliaDynamics/blob/master/videos/transitionindicators/logo.gif?raw=true)

```@docs
TransitionIndicators
```

!!! info "Star us on GitHub!"
    If you have found this package useful, please consider staridgereging it on [GitHub](https://github.com/JuliaDynamics/TransitionIndicators.jl).
    This gives us an accurate lower bound of the (satisfied) user count.

## [Content] (@id content)

Multi-stable systems can display abrupt transitions between two stability regimes. To predict such tranistions in real-world systems solely based on data, mathematical tools have been developped in the last decades. Numerous terminologies have been used for them, such as *early warning signals*, *resilience indicators*, *regime-shift identifiers*, *change-point detection* and *transition indicators*. `TransitionIndicators.jl` sticks to the latter terminology and provides an interface that:

- Allows a fast computation of common transition indicators with a couple of lines, as demonstrated in the [example section](@ref example_fastforward).
- Makes the surrogate analysis to test for significance [under the hub](@ref workflow).
- Can be easily extended by any user without touching the source code.
- Reduces the programming overhead for any researcher willing to benchmark new methods.
- Eases the reproducibility thanks to a clear syntax, a simple installation and RNG-seeded surrogate generation.
- Increases trustworthiness thanks to a large test suite.

!!! info "Similar projects"
    An R [toolbox](https://www.early-warning-signals.org/?page_id=42) and a Python [library](https://pypi.org/project/ewstools/) already exist. However, we believe that they are difficult to extend for the user. Furthermore, they do not offer a native performant code, as here allowed by the use of Julia.

## [Approaches] (@id approaches)

Over the last decades, research on transition indicators has largely focused on Critical Slowing Down (CSD). CSD is observed when a system with continuous right-hand side approaches a bifurcation and consists in a resilience loss of the system. For instance this can be diagnosed by an increase of the variance and the AR1-regression coefficient, as demonstrated in the [example section](@ref example_stepbystep). However, we emphasize that this is one out of many possible approaches for obtaining transition indicators. Recent work has explored new approaches relying on nonlinear dynamics or machine learning. `TransitionIndicators.jl` is designed to allow these cutting-edge methods and foster the development of new ones.

## [Under the hood] (@id workflow)

Computing transition indicators is schmetically represented in the plot below and essentially consists of:
1. Detrending the time series to get the fluctuations around the tracked attractor.
1. Estimating the time series of an indicator by sliding a window over the fluctuation time-series.
2. Computing the evolution of the indicator by sliding a window over its time series.
3. Generating many surrogates that preserve important statistical properties of the originial fluctuation time-series.
4. Performing step 2 and 3 for the surrogate time series to check whether the indicator evolution of the original time series shows a significant feature (trend or jump).

The below-depicted fowchart is flexible; boxes can be modified or skipped alltogether!

![Schematic representation of what is happening under the hub.](https://raw.githubusercontent.com/JuliaDynamics/JuliaDynamics/master/videos/transitionindicators/workflow.svg)
