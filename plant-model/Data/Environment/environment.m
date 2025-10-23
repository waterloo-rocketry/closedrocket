%% Winds

%%% Constant wind
wind_const_direction = deg2rad(180);
wind_const_strength = 12; % m/s

%%% discrete gusts
wind_gust1_start = 5; % s
wind_gust1_length = [240 230 220]; % m
wind_gust1_amplitude = 40; % m/s
wind_gust1_distribution = rand(3,1) .* [0.01; 1; 1]; % factor of gust on each axis

wind_gust2_start = 10;
wind_gust2_length = [400 300 500]; % m
wind_gust2_amplitude = 70; % m/s
wind_gust2_distribution = rand(3,1) .* [0.01; 1; 1]; % factor of gust on each axis

wind_gust3_start = 80;
wind_gust3_length = [400 300 500]; % m
wind_gust3_amplitude = 30; % m/s
wind_gust3_distribution = -rand(3,1) .* [0.01; 1; 1]; % factor of gust on each axis