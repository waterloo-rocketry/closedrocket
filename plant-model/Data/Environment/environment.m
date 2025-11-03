%% Winds

%%% Adding the data set to be available 
%%%T = readtable('plant-model\Data\historical_wind_data\winddata.csv'); 

%get the folder script is in
thisFilePath = mfilename('plant-model\Data\historical_wind_data\winddata.csv');
[thisFolder, ~, ~] = fileparts(thisFilePath);

%build a path to CSV relative to this script
dataPath = fullfile(thisFolder, '..', 'historical_wind_data', 'winddata.csv');

%read the CSV
T = readtable(dataPath, 'VariableNamingRule', 'preserve');


%hand calc stdev here


% Extract relevant wind data from the table
windData = T{:, {'temperature', 'pressure'}};


%%% Constant wind
%wind_const_direction = deg2rad(180);
%wind_const_strength = 10; % m/s

%%% discrete gusts
wind_gust1_start = 10; % s
wind_gust1_length = [120 120 80]; % m
wind_gust1_amplitude = 5; % m/s
wind_gust1_distribution = rand(3,1); % factor of gust on each axis

wind_gust2_start = 20;
wind_gust2_length = [40 30 20]; % m
wind_gust2_amplitude = 10; % m/s
wind_gust2_distribution = rand(3,1); % factor of gust on each axis