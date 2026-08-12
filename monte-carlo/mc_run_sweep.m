function summary = mc_run_sweep(opts)
% Run a Polaris Monte Carlo batch.

    if nargin < 1 || isempty(opts)
        opts = mc_default_options();
    end

    valid_wind_modes = ["none", "discrete", "historical", "mixed"];
    if ~ismember(string(opts.wind.mode), valid_wind_modes)
        error("mc_run_sweep:InvalidWindMode", ...
            "opts.wind.mode must be none, discrete, historical, or mixed.");
    end

    root = char(opts.project_root);
    old_folder = pwd;
    cleanup = onCleanup(@() cd(old_folder));
    cd(root);

    %% Configure the nominal plant
    batch_dir = fullfile(root, char(opts.output_root), ...
        "batch" + string(opts.batch_name));
    if ~exist(batch_dir, "dir")
        mkdir(batch_dir);
    end

    baseline_file = fullfile(batch_dir, "plant_model_baseline.mat");
    rng(opts.seed + 100000, "twister");
    baseline = save_sim_parameters(baseline_file);

    %% Create the sampled simulation inputs
    number_simulations = opts.number_simulations;
    samples = repmat(struct( ...
        "vars", struct(), "meta", struct(), "parameters", struct()), ...
        1, number_simulations);
    simin(1:number_simulations) = ...
        Simulink.SimulationInput(opts.model_name);

    for i = 1:number_simulations
        rng(opts.seed + i - 1, "twister");
        simin(i) = setModelParameter(simin(i), ...
            "StopTime", num2str(opts.stop_time), "CaptureErrors", "on");
        simin(i) = simin(i).loadVariablesFromMATFile(baseline_file);

        samples(i) = mc_sample_inputs(opts, baseline, i);
        sample_names = fieldnames(samples(i).vars);
        for j = 1:numel(sample_names)
            name = sample_names{j};
            simin(i) = simin(i).setVariable( ...
                name, samples(i).vars.(name));
        end
    end

    %% Run the simulations
    close_system(opts.model_name, 0);
    if opts.run_parallel
        simout = parsim(simin, "ShowProgress", char(opts.show_progress));
    else
        for i = 1:number_simulations
            simout(i) = sim(simin(i), ...
                "ShowProgress", char(opts.show_progress)); %#ok<AGROW>
        end
    end
    close_system(opts.model_name, 0);

    %% Save each run and identify errors or unstable results
    error_id = [];
    error_messages = strings(1, number_simulations);
    unstable_id = [];

    for k = 1:number_simulations
        filename = fullfile(batch_dir, sprintf("sim_%d.mat", k));
        parameter_overrides = struct( ...
            "sampled", samples(k).parameters, ...
            "modelVariables", samples(k).vars);

        if ~isempty(simout(k).ErrorMessage)
            error_id(end + 1) = k; %#ok<AGROW>
            error_messages(k) = string(simout(k).ErrorMessage);
            save_failed_run(simout(k), filename, parameter_overrides, ...
                baseline_file, opts);
            continue;
        end

        log = save_log(simout(k), filename, opts.sample_rate_hz, ...
            "ParameterOverrides", parameter_overrides, ...
            "ParameterFile", baseline_file, ...
            "Mode", "sil", ...
            "ModelName", opts.model_name);

        if isfield(log.signals, "cov_norm") && ...
                any(log.signals.cov_norm(:, end) > opts.P_threshold, "all")
            unstable_id(end + 1) = k; %#ok<AGROW>
        end
    end

    %% Save the batch summary
    error_count = length(error_id);
    error_ratio = error_count / number_simulations;
    unstable_count = length(unstable_id);
    unstable_ratio = unstable_count / number_simulations;

    summary = struct();
    summary.batch_dir = batch_dir;
    summary.baseline_file = baseline_file;
    summary.number_simulations = number_simulations;
    summary.error_id = error_id;
    summary.error_messages = error_messages;
    summary.error_count = error_count;
    summary.error_ratio = error_ratio;
    summary.unstable_id = unstable_id;
    summary.unstable_count = unstable_count;
    summary.unstable_ratio = unstable_ratio;

    filename = fullfile(batch_dir, "result_summary.mat");
    save(filename, ...
        "number_simulations", "error_id", "error_messages", ...
        "error_count", "error_ratio", ...
        "unstable_id", "unstable_count", "unstable_ratio", ...
        "samples", "opts", "summary");

    if opts.plot.enable
        mc_plot_batch(opts);
    end
end

function save_failed_run(simout, filename, parameter_overrides, ...
        baseline_file, opts)
% Save diagnostics and any partial signals produced by a failed run.

    log = [];
    logging_error = "";

    try
        log = save_log(simout.logsout, filename, opts.sample_rate_hz, ...
            "ParameterOverrides", parameter_overrides, ...
            "ParameterFile", baseline_file, ...
            "Mode", "sil", ...
            "ModelName", opts.model_name);
    catch err
        logging_error = string(err.message);
    end

    if isempty(log)
        [~, name, extension] = fileparts(filename);
        relative_folder = fullfile(char(opts.output_root), ...
            "batch" + string(opts.batch_name));

        log = struct();
        log.format = "closedrocket-log";
        log.formatVersion = 1;
        log.createdOn = datetime("now");
        log.filePath = string(relative_folder);
        log.fileName = string(name) + string(extension);
        log.sampleRateHz = opts.sample_rate_hz;
        log.stopTime = NaN;
        log.time = zeros(0, 1);
        log.signalNames = strings(0, 1);
        log.signals = struct();
        log.mode = "sil";
        log.modelName = string(opts.model_name);
        log.parameterFile = fullfile(relative_folder, ...
            "plant_model_baseline.mat");
        log.parameterOverrides = parameter_overrides;
    end

    log.status = "failed";
    log.errorMessage = string(simout.ErrorMessage);
    log.partialSignalsAvailable = ~isempty(fieldnames(log.signals));
    if strlength(logging_error) > 0
        log.loggingError = logging_error;
    end

    simulation_output = simout;
    try
        save(filename, "log", "simulation_output", "-v7.3");
    catch err
        warning("mc_run_sweep:SimulationOutputNotSaved", ...
            "Could not save the raw output for failed run '%s': %s", ...
            filename, err.message);
        save(filename, "log", "-v7.3");
    end
end
