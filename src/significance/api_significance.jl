"""
    TransitionsSignificance

Supertype used to test for significance in [`significant_transitions`](@ref).
Valid subtypes are:

- [`SurrogatesSignificance`](@ref).
"""
abstract type TransitionsSignificance end


"""
    significant_transitions(res::IndicatorsChangesResults, signif::TransitionsSignificance)

Estimate significant transtions in `res` using the method described by `signif`.
Return `flags`, a Boolean matrix with identical size as `res.x_change`.
It contains trues wherever a change metric of `res` is deemed significant.
"""
function significant_transitions(::IndicatorsChangesResults, ::TransitionsSignificance) end