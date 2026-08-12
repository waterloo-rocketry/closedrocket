function opts = mc_default_options()
%MC_DEFAULT_OPTIONS Default settings for Polaris Monte Carlo batches.

    opts = struct();
    opts.project_root = fileparts(fileparts(mfilename("fullpath")));
    opts.model_name = "plant-model/CC_Flight_Simulation";
    opts.batch_name = "_polaris_mc";
    opts.output_root = "monte-carlo";
    opts.number_simulations = 100;
    opts.stop_time = 220;
    opts.P_threshold = 5000;
    opts.seed = 42026;
    opts.run_parallel = true;
    opts.show_progress = "on";
    opts.sample_rate_hz = 100;

    opts.plot = struct();
    opts.plot.enable = true;
    opts.plot.export = true;
    opts.plot.visible = "on";
    opts.plot.percentiles = [80, 95];
    opts.plot.time_limits = [0, 80];
    opts.plot.max_runs = inf;

    opts.wind = struct();
    opts.wind.mode = "mixed"; % "none", "discrete", "historical", or "mixed"
    opts.wind.historical_probability = 0.5;
    opts.wind.layer_threshold_mps = 100;
    opts.wind.const_strength_mps = struct("mean", 10, "range", [0, 20]);
    opts.wind.gust_amplitude_mps = struct("mean", 10, "range", [0, 20]);

    opts.geometry = struct();
    opts.geometry.canard_cant_zero_deg = struct( ...
        "mean", 0.2, "range", [-2, 2]);
    opts.geometry.canard_spine_effectiveness = struct( ...
        "mean", 0.5, "range", [0.3, 1]);
    opts.geometry.cp_shift_m = struct( ...
        "mean", 0, "range", [-0.1, 0.1]);
    opts.geometry.canard_radial_cp_offset_m = struct( ...
        "mean", 0, "range", [-0.02, 0.02]);

    opts.performance = struct();
    opts.performance.engine_thrust_factor = struct( ...
        "mean", 1, "range", [0.8, 1.1]);
    opts.performance.act_backlash_deg = struct( ...
        "mean", 0.2, "range", [0, 0.6]);

    opts.roll_reversal = struct();
    opts.roll_reversal.probability = 0.5;
    opts.roll_reversal.factor_range = [-3, 1.5];
    opts.roll_reversal.time_after_takeoff_s = [0, 40];

end
