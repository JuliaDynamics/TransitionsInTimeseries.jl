"""
    ChangesConfig

Supertype for how "changes" in a timeseries are estimated in [`estimate_changes`](@ref).
"Changes" deemed statistically significant in [`significant_transitions`](@ref)
are "transitions" in the timeseries.

Existing subtypes of `ChangesConfig` are:

 - [`SlidingWindowConfig`](@ref).
 - [`SegmentedWindowConfig`](@ref).
 - [`SlopeChangeConfig`](@ref).
"""
abstract type ChangesConfig end

"""
    estimate_changes(config::ChangesConfig, x [,t]) → result

Estimate possible transitions for input timeseries `x` using the approach specified
in the configuration type `config`, see [`ChangesConfig`](@ref) for possibilities.
`t` is the time vector corresponding to `x`, which defaults to `eachindex(x)`.

Return the output as subtype of [`ChangesResults`](@ref).
The particular form of the output depends on the `config` and is described in its docstring.
Regardless of type, `result` can always be given to
[`significant_transitions`](@ref) to deduce which possible transitions are statistically
significant using a variety of significance tests.
"""
function estimate_changes end
# The function is extended via multiple dispatch in the specific files

"""
    ChangesResults

Supertype used to gather results of [`estimate_changes`](@ref).
The concrete subtype instances are described in the docstrings of configuration types.
"""
abstract type ChangesResults end