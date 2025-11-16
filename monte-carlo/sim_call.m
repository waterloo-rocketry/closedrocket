%% Configure
clear 
folderpath = 'monte-carlo/single_controlled/';
run('configure_plant_model');
save(append(folderpath, 'plant_model_baseline.mat'));
clearvars -except folderpath

model_name = 'plant-model/CC_Flight_Simulation';
simin = Simulink.SimulationInput(model_name);

%%% Stop time
% 55 is apogee, 240 is after main deploy
simin = setModelParameter(simin,"StopTime","100");

%%% Load parameters
simin = simin.loadVariablesFromMATFile(append(folderpath, 'plant_model_baseline.mat'));


%% Run Sim

% clear(get_param('CC_Flight_Simulation','ModelWorkspace'))
simout = sim(simin, 'ShowProgress', 'on');

%% Post processing
[sdt, sdt_vars] = sim_postprocessor(simout);
[in_vars] = sim_postprocessor_in(simin, load(append(folderpath, 'plant_model_baseline.mat')));

%% Save 
save(append(folderpath, 'result.mat'), "sdt", "sdt_vars", 'in_vars');
% load(append(folderpath, 'result.mat'), "sdt", "sdt_vars");


%% Plots

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')


% plot_animation(sdt_vars);

f_sim = figure(1);
plot_state(sdt.rocket_dt, "", 'off');

f_est = figure(2);
plots_1 = plot_state(sdt.est, "", 'off');
% plot_state(sdt.rocket_dt, "\_sim", 'off', plots_1);

f_err = figure(3);
plot_state(sdt.error, "");
% sgtitle("Estimation Error")

% exportgraphics(f_sim, append(folderpath, 'state_sim.pdf'));
% exportgraphics(f_est, append(folderpath, 'state_est.pdf'))
% exportgraphics(f_err, append(folderpath, 'state_err.pdf'))

f_cmd = figure(4);
stairs(sdt.control.Time, rad2deg([sdt.control.(3)]))
% legend("Reference", "Roll angle", "Command")
% ylabel("Angle [deg]")
title("Command [deg]",'FontWeight','Normal')
grid on
% exportgraphics(f_cmd, append(folderpath, 'ctrl_cmd.pdf'))

% figure(5)
% sdt_array{1} = sdt;
% plot_stats_covariance(sdt_array, 'P_norm', 'Covariance stats', [50 90])

set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')