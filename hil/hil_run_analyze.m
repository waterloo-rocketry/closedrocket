%% Configure
clear 
run('configure_plant_model');
save('results/hil/plant_model_baseline.mat');
clear

model_name = 'hil/hil';
simin = Simulink.SimulationInput(model_name);

%%% Stop time
%%% 55 is apogee, 240 is after main deploy
simin = setModelParameter(simin,"StopTime","120");

%%% Load parameters
simin = simin.loadVariablesFromMATFile('results/hil/plant_model_baseline.mat');
% simin = simin.setVariable('wind_const_strength', 5);
% simin = simin.setVariable('canard_cant_zero', 0.1);
% simin = simin.setVariable('engine_thrust_factor', 1);
% simin = simin.setVariable('canard_roll_reversal_factor', 1);


%% Run Sim

% clear(get_param('CC_Flight_Simulation','ModelWorkspace'))
out = sim(simin, 'ShowProgress', 'on');

%% Post processing

[sdt, sdt_vars] = hil_postprocessor(out);

%% Save 

% save("results/hil/result.mat", "sdt", "in_vars");
save("results/hil/result.mat", "sdt");


%% Plots

load("results/hil/result.mat", "sdt");

% timecode_overwrite = seconds(0:0.01:seconds(sdt.rocket_dt.Time(end)));
% sdt.est = retime(sdt.est, timecode_overwrite, 'linear');
% sdt.rocket_dt = retime(sdt.rocket_dt, timecode_overwrite, 'linear');
% sdt.error = retime(sdt.error, timecode_overwrite, 'linear');

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')

f_sim = figure(1);
plot_state_old(sdt.rocket_dt, [0, 50, 0], 1:6);
% sgtitle("Simulation")


f_est = figure(2);
plots_1 = plot_state_old(sdt.est, [0, 50, 0], 1:6);
% sgtitle("Estimation")


f_err = figure(3);
plot_state_old(sdt.error, [0, 50, 0], 1:6);
% sgtitle("Estimation Error")


fontsize(f_sim, 12, "points")
fontsize(f_est, 12, "points")
fontsize(f_err, 12, "points")

% set(f_sim,'units','centimeters','position',[1,1,25,18])
% set(f_est,'units','centimeters','position',[1,1,25,18])
% set(f_err,'units','centimeters','position',[1,1,25,18])


set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')

exportgraphics(f_sim, 'results/hil/sim_hil_sim.pdf', 'ContentType', 'vector')
exportgraphics(f_est, 'results/hil/sim_hil_est.pdf', 'ContentType', 'vector')
exportgraphics(f_err, 'results/hil/sim_hil_err.pdf', 'ContentType', 'vector')

% figure(4)
% title("Control")
% plot(sdt.control.Time, rad2deg(sdt.control.(1)) )
% legend("Command")
% ylabel("Angle [deg]")


% figure(5)
% sdt_array{1} = sdt;
% plot_stats_covariance(sdt_array, 'P_norm', 'Covariance stats', [50 90])
