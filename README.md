# TransitionsInTimeseries.jl

[![](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://JuliaDynamics.github.io/TransitionsInTimeseries.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaDynamics.github.io/TransitionsInTimeseries.jl/stable)
[![CI](https://github.com/JuliaDynamics/TransitionsInTimeseries.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/TransitionsInTimeseries.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/TransitionsInTimeseries.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/TransitionsInTimeseries.jl)
[![Package Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/TransitionsInTimeseries)](https://pkgs.genieframework.com?packages=TransitionsInTimeseries)

A Julia package for estimating transitions (from one dynamic regime or stable state to another) in timeseries and testing the statistical significance of found transitions. It integrates with the entire Julia ecosystem of timeseries analysis, and hence offers thousands of metrics that can indicate transitions right out of the box. It also offers a variety of analysis pipelines for identifying transitions and a variety of statistical pipelines for testing for significance.

In contrast to other existing software with similar target application, TransitionsInTimeseries.jl defines a generic interface for how to find transitions and how to test for significance. Within this interface, it is easy to expand the software in three orthogonal ways:

1. Add new indicators that work with the existing analysis pipelines for finding transitions
2. Add new analysis pipelines for finding transitions
3. Add new ways for testing whether already found transitions are significant

This package is currently under active development and not yet registered in the Julia general registry. To install it, first go into package-manager mode in the Julia REPL (press `]`) and then run
```
add https://github.com/JuliaDynamics/TransitionsInTimeseries.jl
```

All further information is provided in the documentation, which you can either find [online](https://juliadynamics.github.io/TransitionsInTimeseries.jl/dev/) or build locally by running the `docs/make.jl` file.

_Alternative names for this package could have been: Early Warning Signals / Resilience Indicators / Regime-Shift Identifiers / Change-Point Detectors, or however else you want to call them!_