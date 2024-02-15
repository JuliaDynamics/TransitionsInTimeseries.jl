# Developer's documentation

This documentation addresses users that would like to contribute to the software by either solving bugs, improving documentation, or adding new methods.
All contributions come in the form of Pull Requests, for which we strongly advise to follow [good practices in scientific code](https://github.com/JuliaDynamics/GoodScientificCodeWorkshop).

## New indicators or change metrics

As explained already in e.g., [`SlidingIndicatorConfig`](@ref), new indicators or change metrics are standard Julia functions, so you only need to define such a function.

## New pipeline for estimating changes

This means to contribute a fundamentally new "pipeline" for estimating/detecting
transitions in timeseries. This new pipeline defines what a "transition" means.
To add a new pipeline follow these steps:

