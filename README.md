# TransitionsInTimeseries.jl

[![](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://JuliaDynamics.github.io/TransitionsInTimeseries.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaDynamics.github.io/TransitionsInTimeseries.jl/stable)
[![CI](https://github.com/JuliaDynamics/TransitionsInTimeseries.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/TransitionsInTimeseries.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/TransitionsInTimeseries.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/TransitionsInTimeseries.jl)
[![Package Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/TransitionsInTimeseries)](https://pkgs.genieframework.com?packages=TransitionsInTimeseries)

A Julia package that can estimate indicators of transitions (from one dynamic regime or stable state to another) in timeseries. Also bundles the indicators with significance testing via surrogate analysis using [TimeseriesSurrogates.jl](https://github.com/JuliaDynamics/TimeseriesSurrogates.jl). Alternative names for this package could have been: Early Warning Signals / Resilience Indicators / Regime-Shift Identifiers / Change-Point Detectors, or however else you want to call them!

This package is currently under active development and not yet registered in the Julia general registry. To install it, first go into package-manager mode in the Julia REPL (press `]`) and then run
```
add https://github.com/JuliaDynamics/TransitionsInTimeseries.jl
```

All further information is provided in the documentation, which you can either find [online](https://juliadynamics.github.io/TransitionsInTimeseries.jl/dev/) or build locally by running the `docs/make.jl` file.
