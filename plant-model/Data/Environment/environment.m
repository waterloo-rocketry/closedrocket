%% Winds

%%% Constant wind
wind_const_direction = deg2rad(360*rand(1));
wind_const_strength = 12; % m/s

%%% discrete gusts
wind_gust1_start = 5; % s
wind_gust1_length = [240 230 220]; % m
wind_gust1_amplitude = 30; % m/s
wind_gust1_distribution = (2*rand(3,1)-1) .* [0.01; 0.9; 0.5]; % factor of gust on each axis

wind_gust2_start = 15;
wind_gust2_length = [400 300 500]; % m
wind_gust2_amplitude = 30; % m/s
wind_gust2_distribution = (2*rand(3,1)-1) .* [-0.01; -0.5; -0.8]; % factor of gust on each axis

wind_gust3_start = 80;
wind_gust3_length = [400 300 500]; % m
wind_gust3_amplitude = 10; % m/s
wind_gust3_distribution = (2*rand(3,1)-1) .* [0.01; 1; 1]; % factor of gust on each axis