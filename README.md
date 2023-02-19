# TransitionIndicators.jl

[![](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://JuliaDynamics.github.io/TransitionIndicators.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaDynamics.github.io/TransitionIndicators.jl/stable)
[![CI](https://github.com/JuliaDynamics/TransitionIndicators.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/TransitionIndicators.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/TransitionIndicators.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/TransitionIndicators.jl)
[![Package Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/TransitionIndicators)](https://pkgs.genieframework.com?packages=TransitionIndicators)

A Julia package that can estimates indicators of transitions (from one dynamic regime or stable state to another) in timeseries. Also bundles the indicators with significance testing via surrogate analysis using [TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl). Alternative names for this package could have been: Early Warning Signals / Resilience Indicators / Regime Shift Identifiers / Change Point Detectors, or however else you want to call them!

To install it, run `import Pkg; Pkg.add("TransitionIndicators")`.

All further information is provided in the documentation, which you can either find [online](https://juliadynamics.github.io/TransitionIndicators.jl/dev/) or build locally by running the `docs/make.jl` file.
