clear

%% Launch site
elevation = 420; % m above sea level

%% Rocket body
m = 54.9; %mass in kg
Jx = 0.46; % inertia roll
Jy = 49.5; % inertia pitch, yaw
J = diag([Jx, Jy, Jy]);
Jinv = inv(J);

length_cg = 0; % center of gravity
length_cp = -0.5; % center of pressure
area_reference = pi*(0.203/2)^2; % cross section of body tube
c_aero = area_reference * (length_cp-length_cg);
Cn_alpha = 10; % pitch forcing coefficent 
Cn_omega = 0; % pitch damping coefficent 

%% Canards
area_canard = 2 * 0.102 * 0.0508 / 2; % total canard area 
length_canard = 0.203/2 + 0.0508/3; % lever arm of canard to x-axis 
c_canard = area_canard*length_canard; % moment arm * area of canard

%% Environment
g = [-9.81; 0; 0]; % gravitational acceleration in the geographic inertial frame

%% Sensors
%%% S: rotation transform from sensor frame to body frame
%%% d: center of sensor frame relative to body frame

S_board = [1 0, 0;
           0, 1, 0;
           0, 0, 1]; % Onboard STM IMU
S_mti = [1, 0, 0;
         0, 1, 0;
         0, 0, 1]; % Movella MTi
S_ad = [1, 0, 0;
         0, 1, 0;
         0, 0, 1]; % AD breakout board

d_board = [0; 0; 0]; % Onboard STM IMU
d_mti = [0; 0; 0]; % Movella MTi
d_ad = [0; 0; 0]; % AD breakout board

%% save and export
save("gnc/model_params.mat");