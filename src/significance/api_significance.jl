"""
    SignificanceConfig

Supertype used to test for significance in [`significant_transitions`](@ref).
Changes that are statistically significant are "transitions".

Valid subtypes are:

- [`SurrogatesSignificance`](@ref).
- [`SigmaSignificance`](@ref).
- [`QuantileSignificance`](@ref).
- [`ThresholdSignificance`](@ref).
"""
abstract type SignificanceConfig end


"""
    significant_transitions(res::ChangesResults, signif::SignificanceConfig)

Estimate significant transtions in `res` using the method described by `signif`.
Return `flags`, a Boolean matrix with identical size as the changes
stored in `res` (which typically is stored in the field `res.x_change`).
`flags` is `true` wherever a change metric of `res` is deemed significant.
"""
function significant_transitions(::ChangesResults, ::SignificanceConfig) end