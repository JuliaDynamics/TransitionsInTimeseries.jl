module TransitionVisualizations

using TransitionsInTimeseries, Makie

# Default options for plotting utilities
const default_accent_linewidth = 3
const default_colors = ["#7143E0", "#0A9A84", "#191E44", "#AF9327", "#701B80", "#2E6137"]

default_indicator_label(res::IndicatorsChangesResults) = [string(nameof(
    res.config.indicators[i])) for i in eachindex(res.config.indicators)]

function default_chametric_label(res::IndicatorsChangesResults)
    labels = String[]
    for i in eachindex(res.config.change_metrics)
        if res.config.change_metrics[i] isa PrecomputedRidgeRegressionSlope ||
            res.config.change_metrics[i] isa RidgeRegressionSlope
            push!(labels, "Regression slope")
        else
            push!(labels, string(nameof(res.config.change_metrics[i])))
        end
    end
    return labels
end

"""
    TransitionVisualization
Contains:
 - `fig::Figure` the figure,
 - `laxs::Vector{Axis}` the left axes,
 - `raxs::Vector{Axis}` the right axes,
 - `res::IndicatorsChangesResults` the results of the input signal,
 - `colors::Vector` a vector containing the colors to use for plotting.
"""
struct TransitionVisualization{R<:IndicatorsChangesResults}
    fig::Makie.Figure
    laxs::Vector{Makie.Axis}
    raxs::Vector{Makie.Axis}
    res::R
    colors::Vector
end

export TransitionVisualization

# Plot results from full analysis
function TransitionsInTimeseries.plot_changes_significance(res, signif; kwargs...)
    tv = plot_indicator_changes(res; kwargs...)
    plot_significance!(tv, signif; kwargs...)
    return tv
end

# Plot results of original input signal
function TransitionsInTimeseries.plot_indicator_changes(res::SlidingWindowResults;
        colors = default_colors,
        indicator_names = default_indicator_label(res),
        chametric_names = default_chametric_label(res),
        additional_timeseries = nothing,
        accent_linewidth = default_accent_linewidth,
        kwargs...,
    )

    fig, laxs, raxs, n, config = init_rowwise_visualisation(res, colors, indicator_names,
        chametric_names, additional_timeseries, accent_linewidth)
    lineplot_metrics!(raxs, laxs, n, config, res.t_indicator, res.x_indicator,
        res.t_change, res.x_change, colors, accent_linewidth)

    return TransitionVisualization(fig, laxs, raxs, res, colors)
end

function TransitionsInTimeseries.plot_indicator_changes(res::SegmentedWindowResults;
    colors = default_colors,
    indicator_names = default_indicator_label(res),
    chametric_names = default_chametric_label(res),
    additional_timeseries = nothing,
    accent_linewidth = default_accent_linewidth,
    kwargs...,
    )

    fig, laxs, raxs, n, config = init_rowwise_visualisation(res, colors, indicator_names,
        chametric_names, additional_timeseries, accent_linewidth)
    for k in eachindex(res.t_indicator)
        Makie.lines!(laxs[1], res.t[k], res.x[k], color = colors[1],
            linewidth = accent_linewidth)
        lineplot_metrics!(raxs, laxs, n, config, res.t_indicator[k], res.x_indicator[k],
            res.t_change[k], res.x_change[k, :], colors, accent_linewidth)
    end

    return TransitionVisualization(fig, laxs, raxs, res, colors)
end

# utils for plot_indicator_changes
function init_rowwise_visualisation(res, colors, indicator_names,
    chametric_names, additional_timeseries, accent_linewidth)

    config = res.config
    fig = Makie.Figure(size = (700, 450), fontsize = 12)
    rlabels = vcat([""], indicator_names)
    llabels = vcat(["Input"], chametric_names)
    n = length(llabels)

    rowaspect = 5
    raxs = [Makie.Axis(fig[(i-1)*rowaspect+1:i*rowaspect, 1], ylabel = rlabels[i],
        xticklabelsvisible = false, yaxisposition = :right, ygridvisible = false,
        ylabelcolor = colors[2], yticklabelcolor = colors[2]) for i in 1:n]
    laxs = [Makie.Axis(fig[(i-1)*rowaspect+1:i*rowaspect, 1], ylabel = llabels[i],
        xticklabelsvisible = false, ygridvisible = false,
        ylabelcolor = colors[1], yticklabelcolor = colors[1]) for i in 1:n]
    Makie.linkxaxes!(laxs..., raxs...)

    hidedecorations!(raxs[1])
    Makie.lines!(laxs[1], res.t, res.x, color = colors[1], linewidth = accent_linewidth)

    raxs[end].xticklabelsvisible = true
    raxs[end].xlabel = "Time"
    Makie.rowgap!(fig.layout, 10)

    elements = [LineElement(color = (colors[1], transparency), linewidth = lw) for 
        (lw, transparency) in [(accent_linewidth, 1), (1, 0.5)]]
    labels = ["original signal", "surro signals"]
    width = 0.5
    if length(res.t_indicator[1]) > 1
        elements = vcat(elements, [MarkerElement(marker = :circle, color = colors[1],
            strokecolor = :transparent, markersize = ms) for ms in [10, 5]])
        labels = vcat(labels, ["original change metric", "surro change metrics"])
        width = 1
    end
    Legend(fig[0, 1], elements, labels, nbanks = 4, rowgap = 0, colgap = 10,
        width = Relative(width))

    return fig, laxs, raxs, n, config
end

function lineplot_metrics!(raxs, laxs, n, config, t_ind, x_ind, t_cha, x_cha,
    colors, accent_linewidth)
    for i in 2:n
        j = i-1
        if !isnothing(config.indicators)
            Makie.lines!(raxs[i], t_ind, x_ind[:, j], color = colors[2],
                linewidth = accent_linewidth)
        end
        if length(t_cha) > 1
            lines!(laxs[i], t_cha, x_cha[:, j], color = colors[1],
                linewidth = accent_linewidth)
        else
            Makie.scatter!(laxs[i], t_cha, x_cha[j], color = colors[1],
                markersize = 10)
        end
    end
end

# Plot results of surrogates
function TransitionsInTimeseries.plot_significance!(
    tv::TransitionVisualization{<:SlidingWindowResults},
    signif::SurrogatesSignificance;
    flags = nothing,
    nsurro = 20,
    kwargs...,
    )

    (; fig, laxs, raxs, res, colors) = tv
    config = res.config
    lines_over_surro!(raxs, laxs, nsurro, res.t, res.t_indicator, res.t_change, res.x,
        signif, config, flags, colors)
    return tv
end

function TransitionsInTimeseries.plot_significance!(
    tv::TransitionVisualization{<:SegmentedWindowResults},
    signif::SurrogatesSignificance;
    flags = nothing,
    nsurro = 20,
    kwargs...,
    )

    (; fig, laxs, raxs, res, colors) = tv
    config = res.config
    for k in eachindex(res.t_indicator)
        if isnothing(flags)
            lines_over_surro!(raxs, laxs, nsurro, res.t[res.i1[k]:res.i2[k]],
                res.t_indicator[k], res.t_change[k], res.x[res.i1[k]:res.i2[k]],
                signif, config, nothing, colors)
        else
            lines_over_surro!(raxs, laxs, nsurro, res.t[res.i1[k]:res.i2[k]],
                res.t_indicator[k], res.t_change[k], res.x[res.i1[k]:res.i2[k]],
                signif, config, flags[k, :], colors)
        end
    end
    return tv
end

# utils for plot_significance!
function lines_over_surro!(raxs, laxs, nsurro, t, t_ind, t_cha, x, signif, config,
    flags, colors)
    c = zeros(length(config.change_metrics), nsurro)
    for ns in 1:nsurro
        s = TimeseriesSurrogates.surrogate(x, signif.surrogate)
        Makie.lines!(laxs[1], t, s; color = (colors[1], 2/nsurro), linewidth = 1)
        for (i, cha) in enumerate(config.change_metrics)

            if isnothing(config.indicators)
                p = s
            else
                p = windowmap(config.indicators[i], s; width = config.width_ind,
                    stride = config.stride_ind)
                Makie.lines!(raxs[i+1], t_ind, p; color = (colors[2], 2/nsurro),
                    linewidth = 1)
            end
            if length(t_cha) > 1
                q = windowmap(cha, p; width = config.width_cha, stride = config.stride_cha)
                Makie.lines!(laxs[i+1], t_cha, q; color = (colors[1], 2/nsurro),
                    linewidth = 1)
            else
                cha = precompute(cha, eachindex(p))
                q = windowmap(cha, p; width = length(p), stride = 1)
                Makie.scatter!(laxs[i+1], t_cha, q[1], color = (colors[1], 2/nsurro),
                    markersize = 5)
            end
            if !isnothing(flags) && ns == 1 && length(t_cha) > 1
                Makie.vlines!(laxs[i+1], t_cha[flags[:, i]];
                    color = (:black, 0.1))
                elem = [PolyElement(color = (:black, 0.5), strokecolor = :transparent)]
                axislegend(laxs[i+1], elem, ["p<$(signif.p)"], position = :lt)
            elseif length(t_cha) == 1
                c[i, ns] = q[1]
            end
        end
    end
end

end