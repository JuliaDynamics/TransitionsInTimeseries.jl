"""
    IndicatorsChangesConfig

Supertype used to define how indicators and their changes are estimated in
[`estimate_indicator_changes`](@ref). Valid subtypes are:

 - [`SlidingWindowConfig`](@ref).
 - [`SegmentedWindowConfig`](@ref).
"""
abstract type IndicatorsChangesConfig end