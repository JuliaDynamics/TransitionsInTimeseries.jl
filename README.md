# TransitionsInTimeseries.jl

[![](https://img.shields.io/badge/docs-dev-lightblue.svg)](https://JuliaDynamics.github.io/TransitionsInTimeseries.jl/dev)
[![](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaDynamics.github.io/TransitionsInTimeseries.jl/stable)
[![CI](https://github.com/JuliaDynamics/TransitionsInTimeseries.jl/workflows/CI/badge.svg)](https://github.com/JuliaDynamics/TransitionsInTimeseries.jl/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/JuliaDynamics/TransitionsInTimeseries.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JuliaDynamics/TransitionsInTimeseries.jl)
[![Package Downloads](https://shields.io/endpoint?url=https://pkgs.genieframework.com/api/v1/badge/TransitionsInTimeseries)](https://pkgs.genieframework.com?packages=TransitionsInTimeseries)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.06464/status.svg)](https://doi.org/10.21105/joss.06464)

TransitionsInTimeseries.jl is a free and open-source software to easily analyse transitions within timeseries in a reproducible, performant, extensible and reliable way. In contrast to other existing software with similar target application, TransitionsInTimeseries.jl defines a generic interface for how to find transitions and how to test for significance. Within this interface, it is easy to expand the software in three orthogonal ways:

1. Provide the analysis pipelines with new indicators, which can be either self-written or imported from other packages. In particular, the latter offers thousands of metrics that can indicate transitions right out of the box.
2. Add new analysis pipelines for finding transitions.
3. Add new ways for significance testing.

TransitionsInTimeseries is a registered Julia package and can be installed by running:
```
] add TransitionsInTimeseries
```

All further information is provided in the documentation, which you can either find [online](https://juliadynamics.github.io/TransitionsInTimeseries.jl/dev/) or build locally by running the `docs/make.jl` file.

_Alternative names for this package could have been: Early Warning Signals / Resilience Indicators / Regime-Shift Identifiers / Change-Point Detectors, or however else you want to call them!_