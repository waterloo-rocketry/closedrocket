%% Choose rocket
% run('plant-model/rockets/Borealis/borealis.m')
% run('plant-model/rockets/testflight/testflight.m')
% run('plant-model/rockets/Aurora/aurora.m')
run('plant-model/rockets/Polaris/polaris.m')

%%% sensors
run('plant-model/Data/Sensors/sensors_polaris.m')

%%% Settings
chute_enable = [1, 0]; % no recovery is = 0
time_idle = 20; % wait time on the rail before launch

%% environment
% enable wind disturbances
wind_disturbance_enable = 1; % no disturbances is = 0
wind_discrete_enable = 1;   % Toggle for type; historic layered wind is 0, old wind model is = 1
wind_layer_threshold = 100; % threshold windspeed at any altitude where launch would not happen
run('plant-model/Data/Environment/wind_discrete.m')
[wind_heights, wind_vectors] = wind_historic(wind_layer_threshold);

%% data pre-processing
run('plant-model/scripts/data_preparation.m')

