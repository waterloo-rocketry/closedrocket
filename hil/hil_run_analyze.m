%% Configure
clear 
run('configure_plant_model');
% save('hil/results/plant_model_baseline.mat');
% clear

model_name = 'hil/hil';
simin = Simulink.SimulationInput(model_name);

%%% Stop time
%%% 55 is apogee, 240 is after main deploy
% simin = setModelParameter(simin,"StopTime","120");

%%% Load parameters
% simin = simin.loadVariablesFromMATFile('hil/results/plant_model_baseline.mat');
% simin = simin.setVariable('wind_const_strength', 5);
% simin = simin.setVariable('canard_cant_zero', 0.1);
% simin = simin.setVariable('engine_thrust_factor', 1);
% simin = simin.setVariable('canard_roll_reversal_factor', 1);


%% Run Sim

% clear(get_param('CC_Flight_Simulation','ModelWorkspace'))
out = sim(simin, 'ShowProgress', 'on');

%% Post processing

[sdt, sdt_vars] = hil_postprocessor(out);
% [in_vars] = sim_postprocessor_in(simin, load('hil/results/plant_model_baseline.mat'));

%% Save 

% save("hil/results/result.mat", "sdt", "in_vars");
save("hil/results/result.mat", "sdt");


%% Plots

load("hil/results/result.mat", "sdt");

% timecode_overwrite = seconds(0:0.01:seconds(sdt.rocket_dt.Time(end)));
% sdt.est = retime(sdt.est, timecode_overwrite, 'linear');
% sdt.rocket_dt = retime(sdt.rocket_dt, timecode_overwrite, 'linear');
% sdt.error = retime(sdt.error, timecode_overwrite, 'linear');

figure(1)
plot_state(sdt.rocket_dt);
sgtitle("Simulation")

figure(2)
plots_1 = plot_state(sdt.est);
sgtitle("Estimation")

% figure(3)
% plot_state(sdt.error);
% sgtitle("Estimation Error")

% figure(4)
% title("Control")
% stairs(sdt.control.Time, rad2deg([sdt.control.(1), sdt.control.(4), sdt.control.(3)]))
% legend("Reference", "Roll angle", "Command")
% ylabel("Angle [deg]")


% figure(5)
% sdt_array{1} = sdt;
% plot_stats_covariance(sdt_array, 'P_norm', 'Covariance stats', [50 90])