"""
    IndicatorsChangesConfig

Supertype used to define how indicators and their changes are estimated in
[`estimate_indicator_changes`](@ref). Valid subtypes are:

 - [`SlidingWindowConfig`](@ref).
 - [`SegmentedWindowConfig`](@ref).
"""
abstract type IndicatorsChangesConfig end

"""
    estimate_indicator_changes(config::IndicatorsChangesConfig, x [,t]) → result

Estimate possible transitions for input timeseries `x` using the approach specified
in the configuration type `config`, see [`IndicatorsChangesConfig`](@ref) for possibilities.
`t` is the time vector corresponding to `x`, which defaults to `eachindex(x)`.

Return the output as subtype of [`IndicatorsChangesResults`](@ref).
The particular form of the output depends on the `config` and is described in its docstring.
Regardless of type, `result` can always be given to
[`significant_transitions`](@ref) to deduce which possible transitions are statistically
significant using a variety of significance tests.
"""
function estimate_indicator_changes end
# The function is extended via multiple dispatch in the specific files

"""
    IndicatorsChangesResults

Supertype used to gather results of [`estimate_indicator_changes`](@ref).
The concrete subtype instances are described in the docstrings of configuration types.
"""
abstract type IndicatorsChangesResults end