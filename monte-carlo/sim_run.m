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

