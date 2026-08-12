%% Polaris Monte Carlo batch

opts = mc_default_options();


opts.batch_name = "_polaris_mc";
opts.number_simulations = 1000;
opts.stop_time = 120;
opts.P_threshold = 5000;
opts.seed = 42026;

% Wind modes: "none", "discrete", "historical", or "mixed".
opts.wind.mode = "mixed";
opts.wind.historical_probability = 0.5;

summary = mc_run_sweep(opts);
disp(summary);
