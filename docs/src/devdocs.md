# Developer's documentation

This documentation addresses users that would like to contribute to the software by either solving bugs, improving documentation, or adding new methods.
All contributions come in the form of Pull Requests, for which we strongly advise to follow [good practices in scientific code](https://github.com/JuliaDynamics/GoodScientificCodeWorkshop), which means that they properly formatted, documented, tested, etc.

## New indicators or change metrics

As explained already in e.g., [`SlidingIndicatorConfig`](@ref), new indicators or change metrics are standard Julia functions, so you only need to define such a function (and document it, test it, etc.).

## New pipeline for estimating changes

This means to contribute a fundamentally new "pipeline" for estimating/detecting
transitions in timeseries. This new pipeline defines what a "transition" means.
To add a new pipeline follow these steps:

1. Define a new subtype of [`ChangesConfig`](@ref)
2. Define a new subtype of [`ChangesResults`](@ref)
3. Add a method for [`estimate_changes`](@ref) which accepts
   the new `ChangesConfig` subtype you defined and
   returns the `ChangesResults` subtype you defined.

And that's it!

_Optionally_ you can further extend:

- `TransitionsInTimeseries.plot_indicator_changes` with the new `ChangesResults` type you defined. Note the plotting functions are in the `ext` folder of the repository.
- `significant_transitions(res::ChangesResults, signif::SignificanceConfig)`
  with your new `ChangesResults` type and as many `SignificanceConfig`
  subtypes you have the capacity to extend for.


## New pipeline for estimating significance

Statistically significant changes are "transitions".
However, what "significant" means is not universal. There are different ways to
test for significance and TransitionsInTimeseries.jl allows various methods.

To add a new pipeline for estimating significance follow these steps:

1. Define a new subtype of [`SignificanceConfig`](@ref)
2. Extend `significant_transitions(res::ChangesResults, signif::SignificanceConfig)`
   with your new type and as many `ChangesResults` subtypes as you have
   capacity for.

And that's it! Unfortunately so far we have found no way to make the
`significant_transitions` code agnostic of the changes result, so you would
need to add a method manually for every `ChangesResults`.

_Optionally_ you can further extend `TransitionsInTimeseries.plot_significance!`
(which you can find in the `ext` folder) for visualizing the significance results.
