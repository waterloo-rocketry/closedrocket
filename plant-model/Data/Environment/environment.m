%% Winds

%%% Constant wind
% wind_const_direction = deg2rad(180);
% wind_const_strength = 10; % m/s

%%% discrete gusts
wind_gust1_start = 10; % s
wind_gust1_length = [120 120 80]; % m
wind_gust1_amplitude = 5; % m/s
wind_gust1_distribution = rand(3,1); % factor of gust on each axis

wind_gust2_start = 20;
wind_gust2_length = [40 30 20]; % m
wind_gust2_amplitude = 10; % m/s
wind_gust2_distribution = rand(3,1); % factor of gust on each axis


%%% windspeed wrt height from csv
wind_data = readtable("plant-model/Data/environment/sim_parameters_historical_aug2023+aug2024_no-clouds_2025-06-03.csv", ...
    ReadVariableNames=true, VariableNamingRule="preserve");

%% pick a random timestamp and process
row_idx = randi(height(wind_data));
% row_idx = 1;
row = wind_data(row_idx, :);
wind_time = row.date;
fprintf('Selected timestamp: %s\n', wind_time);

cols = wind_data.Properties.VariableNames;

% Speed columns (numeric names, e.g., '5', '10', etc.)
speed_cols = cols(~cellfun(@isempty, regexp(cols, "\d+$")));

% Direction columns (e.g., 'direction [5m]', 'direction [10m]')
dir_cols = cols(~cellfun(@isempty, regexp(cols, "direction")));

% Convert column names to numeric heights [m]
heights = str2double(speed_cols);

% Extract wind speeds and directions for selected timestamp
speeds = row{1, speed_cols}; % [m/s]
dirs   = row{1, dir_cols};   % [deg]

%% Convert to 3D wind vectors
% im assuming rocket axes (x = up, y = north, z = east)
% and that the direction is deg from north
wind_vectors = zeros(length(heights), 3);

for i = 1:length(heights)
    deg = dirs(i);
    rad = deg2rad(deg);
    % stolen from the way wind_const_strength and wind_const_direction were used before
    wind_vectors(i, :) = speeds(i) * [0, cos(rad), sin(rad)]; % no vertical component
end

% N x 1
wind_heights = heights';
% N x 3
wind_vectors = wind_vectors;


% fprintf('heights:\n');
% disp(wind_heights);
% fprintf("windvectors:\n");
% disp(wind_vectors)