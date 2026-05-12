clear

%% Launch site
elevation = 420; % [m], launch altitude above sea level

%% Rocket body
m = 54.9; % [kg], rocket dry mass 
Jx = 0.46; % [kg m^2] rocket dry roll inertia
Jy = 49.5; % [kg m^2] rocket dry pitch, yaw inertia
J = diag([Jx, Jy, Jy]);
Jinv = inv(J); % precompute inverse inertia

length_cg = 0; % [m], center of gravity
length_cp = -0.5; % [m], center of pressure (negative is aft of cg)
area_reference = pi*(0.203/2)^2; % [m^2], cross section of body tube
c_aero = area_reference * (length_cp-length_cg); % precompute aero constants
Cn_alpha = 10; % [-], pitch forcing coefficent 
Cn_omega = 0; % [-], pitch damping coefficent 

%% Canards
area_canard = 2 * 0.102 * 0.0508 / 2; % [m^2], total canard area 
length_canard = 0.203/2 + 0.0508/3; % [m], lever arm of canard to x-axis 
c_canard = area_canard*length_canard; % precompute aero constants

%% Environment
g = [-9.81; 0; 0]; % [m/s^2], gravitational acceleration in the geographic inertial frame

%% Sensors
%%% S: [1], rotation transform from sensor frame to body frame
%%% d: [m], center of sensor frame relative to body frame

S_board = [1 0, 0;
           0, 1, 0;
           0, 0, 1]; % Onboard STM IMU
S_mti = [1, 0, 0;
         0, 1, 0;
         0, 0, 1]; % Movella MTi
S_ad = [1, 0, 0;
         0, 1, 0;
         0, 0, 1]; % AD breakout board

d_board = [1.21; 0; 0]; % Onboard STM IMU
d_mti = [1.21; 0; 0]; % Movella MTi
d_ad = [1.21; 0; 0]; % AD breakout board

%% save and export
save("gnc/model_params.mat");