module TransitionVisualizations

using TransitionsInTimeseries, Makie

# Default options for plotting utilities
const default_accent_linewidth = 3
const default_colors = ["#7143E0", "#0A9A84", "#191E44", "#AF9327", "#701B80", "#2E6137"]

default_indicator_label(res::ChangesResults) = [shortname(
    ind) for ind in res.config.indicators]

default_chametric_label(res::ChangesResults) = [shortname(
    cha_metric) for cha_metric in res.config.change_metrics]

shortname(metric) = string(nameof(metric))
shortname(metric::PrecomputedRidgeRegressionSlope) = "Regression slope"
shortname(metric::RidgeRegressionSlope) = "Regression slope"

# Plot results from full analysis
function TransitionsInTimeseries.plot_changes_significance(res, signif; kwargs...)
    fig = plot_indicator_changes(res; kwargs...)
    plot_significance!(fig, res, signif; kwargs...)
    return fig
end

# Plot results of original input signal
function TransitionsInTimeseries.plot_indicator_changes(res::SlidingWindowResults;
        colors = default_colors,
        indicator_names = default_indicator_label(res),
        chametric_names = default_chametric_label(res),
        accent_linewidth = default_accent_linewidth,
        kwargs...,
    )

    fig, n = init_rowwise_visualisation(res, colors, indicator_names,
        chametric_names, accent_linewidth)
    lineplot_metrics!(fig, n, res.config, res.t_indicator, res.x_indicator,
        res.t_change, res.x_change, colors, accent_linewidth)

    return fig
end

function TransitionsInTimeseries.plot_indicator_changes(res::SegmentedWindowResults;
    colors = default_colors,
    indicator_names = default_indicator_label(res),
    chametric_names = default_chametric_label(res),
    accent_linewidth = default_accent_linewidth,
    kwargs...,
    )

    fig, n = init_rowwise_visualisation(res, colors, indicator_names,
        chametric_names, accent_linewidth)
    for k in eachindex(res.t_indicator)
        Makie.lines!(contents(fig[1, 1])[1], res.t[k], res.x[k], color = colors[1],
            linewidth = accent_linewidth)
        lineplot_metrics!(fig, n, res.config, res.t_indicator[k], res.x_indicator[k],
            res.t_change[k], res.x_change[k, :], colors, accent_linewidth)
    end

    return fig
end

# utils for plot_indicator_changes
function init_rowwise_visualisation(res, colors, indicator_names,
    chametric_names, accent_linewidth)

    fig = Makie.Figure(size = (700, 450), fontsize = 12)
    rlabels = vcat([""], indicator_names)
    llabels = vcat(["Input"], chametric_names)
    n = length(llabels)

    # rowaspect = 5
    raxs = [Makie.Axis(fig[i, 1], ylabel = rlabels[i],
        xticklabelsvisible = false, yaxisposition = :right, ygridvisible = false,
        ylabelcolor = colors[2], yticklabelcolor = colors[2]) for i in 1:n]
    laxs = [Makie.Axis(fig[i, 1], ylabel = llabels[i],
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
    rowsize!(fig.layout, 0, Auto(0.3))
    return fig, n
end

function lineplot_metrics!(fig, n, config, t_ind, x_ind, t_cha, x_cha,
    colors, accent_linewidth)
    for i in 2:n
        j = i-1
        if !isnothing(config.indicators)
            Makie.lines!(contents(fig[i, 1])[2], t_ind, x_ind[:, j], color = colors[2],
                linewidth = accent_linewidth)
        end
        if length(t_cha) > 1
            lines!(contents(fig[i, 1])[1], t_cha, x_cha[:, j], color = colors[1],
                linewidth = accent_linewidth)
        else
            Makie.scatter!(contents(fig[i, 1])[1], t_cha, x_cha[j], color = colors[1],
                markersize = 10)
        end
    end
end

# Plot results of surrogates
function TransitionsInTimeseries.plot_significance!(
    fig::Figure,
    res::SlidingWindowResults,
    signif::SurrogatesSignificance;
    colors = default_colors,
    flags = nothing,
    nsurro = 20,
    kwargs...,
    )

    lines_over_surro!(fig, nsurro, res.t, res.t_indicator, res.t_change, res.x,
        signif, res.config, flags, colors)
    return nothing
end

function TransitionsInTimeseries.plot_significance!(
    fig::Figure,
    res::SegmentedWindowResults,
    signif::SurrogatesSignificance;
    colors = default_colors,
    flags = nothing,
    nsurro = 20,
    kwargs...,
    )

    for k in eachindex(res.t_indicator)
        if isnothing(flags)
            lines_over_surro!(fig, nsurro, res.t[res.i1[k]:res.i2[k]],
                res.t_indicator[k], res.t_change[k], res.x[res.i1[k]:res.i2[k]],
                signif, res.config, nothing, colors)
        else
            lines_over_surro!(fig, nsurro, res.t[res.i1[k]:res.i2[k]],
                res.t_indicator[k], res.t_change[k], res.x[res.i1[k]:res.i2[k]],
                signif, res.config, flags[k, :], colors)
        end
    end
    return nothing
end

# utils for plot_significance!
function lines_over_surro!(fig, nsurro, t, t_ind, t_cha, x, signif, config,
    flags, colors)
    c = zeros(length(config.change_metrics), nsurro)
    for ns in 1:nsurro
        s = TimeseriesSurrogates.surrogate(x, signif.surrogate)
        Makie.lines!(contents(fig[1, 1])[1], t, s; color = (colors[1], 2/nsurro), linewidth = 1)
        for (i, cha) in enumerate(config.change_metrics)

            if isnothing(config.indicators)
                p = s
            else
                p = windowmap(config.indicators[i], s; width = config.width_ind,
                    stride = config.stride_ind)
                Makie.lines!(contents(fig[i+1, 1])[2], t_ind, p; color = (colors[2], 2/nsurro),
                    linewidth = 1)
            end
            if length(t_cha) > 1
                q = windowmap(cha, p; width = config.width_cha, stride = config.stride_cha)
                Makie.lines!(contents(fig[i+1, 1])[1], t_cha, q; color = (colors[1], 2/nsurro),
                    linewidth = 1)
            else
                cha = precompute(cha, eachindex(p))
                q = windowmap(cha, p; width = length(p), stride = 1)
                Makie.scatter!(contents(fig[i+1, 1])[1], t_cha, q[1], color = (colors[1], 2/nsurro),
                    markersize = 5)
            end
            if !isnothing(flags) && ns == 1 && length(t_cha) > 1
                Makie.vlines!(contents(fig[i+1, 1])[1], t_cha[flags[:, i]];
                    color = (:black, 0.1))
                elem = [PolyElement(color = (:black, 0.5), strokecolor = :transparent)]
                axislegend(contents(fig[i+1, 1])[1], elem, ["p<$(signif.p)"], position = :lt)
            elseif length(t_cha) == 1
                c[i, ns] = q[1]
            end
        end
    end
end

end