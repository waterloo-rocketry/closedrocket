%% Configure
batch_name = '_ascent_200';
number_plots = 200;
% exclude = [88, 177]; %indices
% limit_filesize = 4000; %kB
% limit_velocity = 1000;
percentiles = [80, 95];

%% Load statistical
sdt_array = cell(1, number_plots);
filename = sprintf('results/monte-carlo/batch%s/result_summary.mat', batch_name);
load(filename);
%unstable_id
for k = 1:number_plots
    filename = sprintf('results/monte-carlo/batch%s/sim_%d.mat', batch_name, k);
    load(filename);  % loads variables: sdt, in_vars
    % if any(sdt.est.v(:,1) > limit_velocity) 
    if ismember(k, unstable_id) 
        continue    % skip these for now
    end
    sdt_array{k} = sdt;  % store the sdt struct
end


%% Plot statistical

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')


% f_sim = figure(1);
% plot_stats_state(sdt_array, 'rocket_dt', [80, 95], [0, 50], [2, 2]);
% f_est = figure(2);
% plot_stats_state(sdt_array, 'est', 'Estimation', percentiles, [0, 250]);
% f_err = figure(3);
% plot_stats_state(sdt_array, 'error', [80, 95], [0, 50], 1:6);
% f_cov = figure(4);
% plot_stats_covariance(sdt_array, 'P_norm', [80, 95], [0, 50], 3);
f_control = figure(5);
plot_stats_control(sdt_array, 'control', 'Controller', [80, 95]);

% fontsize(f_sim, 12, "points")
% fontsize(f_est, 12, "points")
% fontsize(f_err, 12, "points")
fontsize(f_cov, 12, "points")

set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')

%% save statistical
% set(f_sim, 'Units', 'normalized', 'WindowState', 'maximized');
% exportgraphics(f_sim, sprintf('results/monte-carlo/batch%s/sim_mc_stats_sim%s.pdf', batch_name, batch_name), 'ContentType', 'vector')
% set(f_est, 'Units', 'normalized', 'WindowState', 'maximized');
% exportgraphics(f_est, sprintf('results/monte-carlo/batch%s/sim_mc_stats_est%s.pdf', batch_name, batch_name), 'ContentType', 'vector')
% set(f_err, 'Units', 'normalized', 'WindowState', 'maximized');
% exportgraphics(f_err, sprintf('results/monte-carlo/batch%s/sim_mc_stats_error%s.pdf', batch_name, batch_name), 'ContentType', 'vector')
% set(f_control, 'Units', 'normalized', 'WindowState', 'maximized');
% exportgraphics(f_control, sprintf('results/monte-carlo/batch%s/sim_mc_stats_control%s.pdf', batch_name, batch_name), 'ContentType', 'vector')
% exportgraphics(f_cov, sprintf('results/monte-carlo/batch%s/sim_mc_stats_cov%s.pdf', batch_name, batch_name), 'ContentType', 'vector');
