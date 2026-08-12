function figs = mc_plot_batch(opts)
%MC_PLOT_BATCH Plot statistical envelopes for a Monte Carlo batch.

    if nargin < 1 || isempty(opts)
        opts = mc_default_options();
    end

    batch_dir = fullfile(opts.project_root, char(opts.output_root), ...
        "batch" + string(opts.batch_name));
    summary_data = load(fullfile(batch_dir, "result_summary.mat"));

    number_plots = min( ...
        summary_data.number_simulations, opts.plot.max_runs);
    sdt_array = cell(1, number_plots);

    for k = 1:number_plots
        filename = fullfile(batch_dir, sprintf("sim_%d.mat", k));
        if ~isfile(filename) || ismember(k, summary_data.error_id) || ...
                ismember(k, summary_data.unstable_id)
            continue;
        end

        run_data = load(filename, "log");
        sdt_array{k} = mc_log_to_sdt(run_data.log);
    end

    set(groot, "defaultAxesTickLabelInterpreter", "latex");
    set(groot, "defaultLegendInterpreter", "latex");
    set(groot, "DefaultTextInterpreter", "latex");

    figs.state = figure("Name", "Monte Carlo Rocket State", ...
        "Visible", opts.plot.visible);
    plot_stats_state(sdt_array, "rocket_dt", opts.plot.percentiles, ...
        opts.plot.time_limits, 1:6);

    figs.estimator = figure("Name", "Monte Carlo Estimator Error", ...
        "Visible", opts.plot.visible);
    plot_stats_state(sdt_array, "error", opts.plot.percentiles, ...
        opts.plot.time_limits, 1:4);

    figs.control = figure("Name", "Monte Carlo Control", ...
        "Visible", opts.plot.visible);
    plot_stats_control(sdt_array, "control", "Controller", ...
        opts.plot.percentiles, opts.plot.time_limits);

    if opts.plot.export
        plot_dir = fullfile(batch_dir, "plots");
        if ~exist(plot_dir, "dir")
            mkdir(plot_dir);
        end

        exportgraphics(figs.state, fullfile(plot_dir, "mc_state.pdf"), ...
            "ContentType", "vector");
        exportgraphics(figs.estimator, ...
            fullfile(plot_dir, "mc_estimator_error.pdf"), ...
            "ContentType", "vector");
        exportgraphics(figs.control, fullfile(plot_dir, "mc_control.pdf"), ...
            "ContentType", "vector");
    end
end
