function figs = plot_controller(source, varargin)
%PLOT_CONTROLLER Plot controller dashboard signals from a saved log.

    if nargin < 1 || isempty(source)
        source = "result.mat";
    end

    parser = inputParser;
    parser.KeepUnmatched = true;
    addParameter(parser, "TimeLimits", []);
    addParameter(parser, "SaveFigures", false);
    addParameter(parser, "OutputDir", "");
    parse(parser, varargin{:});

    data = load_log_data(source, "plot_controller");
    set_log_plot_defaults();

    figs = struct();
    figs.controller = plotControllerDashboard(data, parser.Results.TimeLimits);
    figs = remove_empty_figures(figs);

    if parser.Results.SaveFigures
        outputDir = default_output_dir(source, parser.Results.OutputDir);
        export_log_figures(figs, outputDir);
    end
end

function fig = plotControllerDashboard(data, tlim)
    phiTarget = get_log_signal(data, "phi_target");
    phi = get_log_signal(data, "phi");
    wTarget = get_log_signal(data, "w_target");
    w = first_log_column(get_log_signal(data, "w"));
    cmd = get_log_signal(data, "cmd");
    canardAngle = get_log_signal(data, "canard_angle");
    CL = get_log_signal(data, "CL");
    CLEst = get_log_signal(data, "CL_est");

    if all(cellfun(@isempty, {phiTarget, phi, wTarget, w, cmd, canardAngle, CL, CLEst}))
        fig = [];
        return;
    end

    fig = figure("Name", "Controller Dashboard");
    layout = tiledlayout(fig, 2, 2, "TileSpacing", "compact", "Padding", "compact");

    ax = nexttile(layout);
    plot_signal_group(ax, {phiTarget, phi}, ["phi_target", "phi"], ...
        "Roll Angle", "rad", tlim);

    ax = nexttile(layout);
    plot_signal_group(ax, {cmd, canardAngle}, ["cmd", "canard_angle"], ...
        "Canard Angle", "deg", tlim);

    ax = nexttile(layout);
    plot_signal_group(ax, {wTarget, w}, ["w_target", "w"], ...
        "Angular Rates", "rad/s", tlim);

    ax = nexttile(layout);
    plot_signal_group(ax, {CL, CLEst}, ["CL", "CL_est"], ...
        "Roll Control Coefficient derivative", "", tlim);
end
