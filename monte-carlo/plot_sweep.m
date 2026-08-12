%% Plot a Polaris Monte Carlo batch
opts = mc_default_options();
opts.batch_name = "_polaris_mc";
opts.plot.export = true;
opts.plot.visible = "on";
opts.plot.percentiles = [80, 95];
opts.plot.time_limits = [0, 100];

figs = mc_plot_batch(opts);
disp(figs);
