"""
    plot_indicator_changes(res) → (fig, axs)
Return `fig::Figure` and `axs::Matrix{Axis}`, on which `res::ChangesResults`
has been visualised.

## Keyword arguments:
 - `colors = default_colors` sets the colors of the line plots that are to
 be cycled through. The default correspond to the color scheme of Julia Dynamics.
 - `indicator_names = default_indicator_label(res), chametric_names =
 default_chametric_label(res)` sets the labels for the indicators and the change
 metrics, with the default inferring them from the names of `res.config.indicators`
 and `res.config.change_metrics`.
 - `accent_linewidth = 3` sets the line width for the original signals (the
 surrogates have `linewidth = 1`)
"""
function plot_indicator_changes end

"""
    plot_significance!(axs, res, signif)
Update the `axs::Matrix{Axis}` of a figure obtained with `plot_indicator_changes(res)`
with the `signif::SurrogatesSignificance`.

## Keyword arguments:
 - `flags = nothing` provides the significance flags, for instance obtained by
 thresholding the pvalues.
 - `nsurro = 20` sets the number of surrogates to plot.
"""
function plot_significance! end

"""
    plot_changes_significance(res, signif) → (fig, axs)
Return `fig::Figure` and `axs::Matrix{Axis}`, on which `res::ChangesResults`
and `signif::SurrogatesSignificance` have been visualised.
The source code is as simple as:

```julia
fig, axs = plot_indicator_changes(res; kwargs...)
plot_significance!(axs, res, signif; kwargs...)
```

For more information, refer to [`plot_indicator_changes`](@ref) and
[`plot_significance!`](@ref).
"""
function plot_changes_significance end

export plot_indicator_changes, plot_significance!, plot_changes_significance