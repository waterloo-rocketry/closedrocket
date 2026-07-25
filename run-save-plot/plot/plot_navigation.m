function figs = plot_navigation(source, varargin)
%PLOT_NAVIGATION Plot navigation state estimates and covariance.
%
%   plot_navigation("result.mat") opens a saved log created by save_log.
%   plot_navigation(log, "PlotError", true) also plots derived estimation
%   error without requiring error signals to be saved.

    if nargin < 1 || isempty(source)
        source = "result.mat";
    end

    parser = inputParser;
    addParameter(parser, "TimeLimits", []);
    addParameter(parser, "PlotError", false);
    addParameter(parser, "SaveFigures", false);
    addParameter(parser, "OutputDir", "");
    parse(parser, varargin{:});

    data = load_log_data(source, "plot_navigation");
    set_log_plot_defaults();

    tlim = parser.Results.TimeLimits;
    figs = struct();
    figs.navigation = plotNavigationDashboard(data, tlim);
    figs.covariance = plotCovarianceDashboard(data, tlim);
    if parser.Results.PlotError
        figs.error = plotErrorDashboard(data, tlim);
    end
    figs = remove_empty_figures(figs);

    if parser.Results.SaveFigures
        outputDir = default_output_dir(source, parser.Results.OutputDir);
        export_log_figures(figs, outputDir);
    end
end

function fig = plotNavigationDashboard(data, tlim)
    q = get_log_signal(data, "q");
    qEst = get_log_signal(data, "q_est");
    v = get_log_signal(data, "v");
    vEst = get_log_signal(data, "v_est");
    w = get_log_signal(data, "w");
    wEst = get_log_signal(data, "w_est");
    r = get_log_signal(data, "r");
    altEst = get_log_signal(data, "alt_est");

    if all(cellfun(@isempty, {q, qEst, v, vEst, w, wEst, r, altEst}))
        fig = [];
        return;
    end

    fig = figure("Name", "Navigation State Estimates");
    layout = tiledlayout(fig, 2, 2, "TileSpacing", "compact", "Padding", "compact");

    ax = nexttile(layout);
    plot_signal_group(ax, {q, qEst}, ["q", "q\_est"], ...
        "Attitude Quaternion", "", tlim);

    ax = nexttile(layout);
    plot_signal_group(ax, {v, vEst}, ["v", "v\_est"], ...
        "Velocity", "m/s", tlim);

    ax = nexttile(layout);
    plot_signal_group(ax, {w, wEst}, ["w", "w\_est"], ...
        "Angular Rates", "rad/s", tlim);

    ax = nexttile(layout);
    plot_signal_group(ax, {r, altEst}, ["r", "alt\_est"], ...
        "Position", "m", tlim);
end

function fig = plotCovarianceDashboard(data, tlim)
    covNorm = get_log_signal(data, "cov_norm");
    if isempty(covNorm)
        fig = [];
        return;
    end

    fig = figure("Name", "Navigation Covariance");
    ax = axes(fig);
    plot_signal_group(ax, {covNorm}, "cov\_norm", "Covariance Norm", "", tlim);
end

function fig = plotErrorDashboard(data, tlim)
    qErr = signal_difference(get_log_signal(data, "q"), ...
        get_log_signal(data, "q_est"), "q_error");
    vErr = signal_difference(get_log_signal(data, "v"), ...
        get_log_signal(data, "v_est"), "v_error");
    wErr = signal_difference(get_log_signal(data, "w"), ...
        get_log_signal(data, "w_est"), "w_error");
    altErr = signal_difference(first_log_column(get_log_signal(data, "r")), ...
        get_log_signal(data, "alt_est"), "alt_error");

    if all(cellfun(@isempty, {qErr, vErr, wErr, altErr}))
        fig = [];
        return;
    end

    fig = figure("Name", "Navigation Estimation Error");
    layout = tiledlayout(fig, 2, 2, "TileSpacing", "compact", "Padding", "compact");

    ax = nexttile(layout);
    plot_signal_group(ax, {qErr}, "q\_error", "Quaternion Error", "", tlim);

    ax = nexttile(layout);
    plot_signal_group(ax, {vErr}, "v\_error", "Velocity Error", "m/s", tlim);

    ax = nexttile(layout);
    plot_signal_group(ax, {wErr}, "w\_error", "Angular Rate Error", "rad/s", tlim);

    ax = nexttile(layout);
    plot_signal_group(ax, {altErr}, "alt\_error", "Altitude Error", "m", tlim);
end
