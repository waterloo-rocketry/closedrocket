function set_log_plot_defaults()
%SET_LOG_PLOT_DEFAULTS Use LaTeX text rendering for log plots.

    set(groot, "defaultAxesTickLabelInterpreter", "latex");
    set(groot, "defaultLegendInterpreter", "latex");
    set(groot, "DefaultTextInterpreter", "latex");
end
