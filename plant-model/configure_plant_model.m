%% Choose rocket
% run('plant-model/rockets/Borealis/borealis.m')
% run('plant-model/rockets/testflight/testflight.m')
% run('plant-model/rockets/Aurora/aurora.m')
run('plant-model/rockets/Polaris/polaris.m')

%%% Settings
chute_enable = [1, 0]; % no recovery is = 0
time_idle = 10; % wait time on the rail before launch

%%% environment
run('plant-model/Data/Environment/wind_discrete.m')
run('plant-model/Data/Environment/wind_historic.m')
% enable wind disturbances
wind_dist_enable = 1; % no disturbances is = 0
constant_wind_enable = 0;   % Toggle for type; historic layered wind is 0, old wind model is = 1

%%% sensors
run('plant-model/Data/Sensors/sensors_polaris.m')

%% data pre-processing
run('plant-model/scripts/data_preparation.m')

