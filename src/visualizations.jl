function plot_changes_significance(res::IndicatorsChangesConfig, signif::TransitionsSignificance)
    fig = plot_indicator_changes()
    plot_significance!()
    return fig
end

function plot_indicator_changes end

function plot_significance! end