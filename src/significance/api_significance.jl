"""
    TransitionsSignificance

Supertype used to test for significance in [`significant_transitions`](@ref).
Valid subtypes are:

- [`SurrogatesSignificance`](@ref).
- [`SigmaSignificance`](@ref).
- [`QuantileSignificance`](@ref).
- [`ThresholdSignificance`](@ref).
"""
abstract type TransitionsSignificance end


"""
    significant_transitions(res::ChangesResults, signif::TransitionsSignificance)

Estimate significant transtions in `res` using the method described by `signif`.
Return `flags`, a Boolean matrix with identical size as the changes
stored in `res` (which typically is stored in the field `res.x_change`).
`flags` is `true` wherever a change metric of `res` is deemed significant.
"""
function significant_transitions(::ChangesResults, ::TransitionsSignificance) end