clear

%% Launch site
altitude_initial = 420; % [m], launch altitude above sea level

%% Rocket body
m = 66.342; % [kg], rocket dry mass 
Jx = 0.64; % [kg m^2] rocket dry roll inertia
Jy = 189.5; % [kg m^2] rocket dry pitch, yaw inertia
J = diag([Jx, Jy, Jy]);
Jinv = inv(J); % precompute inverse inertia

length_cg = 0; % [m], center of gravity
length_cp = -1.1; % [m], center of pressure (negative is aft of cg)
area_reference = pi*(0.203/2)^2; % [m^2], cross section of body tube
c_aero = area_reference * (length_cp-length_cg); % precompute aero constants
Cn_alpha = 10; % [-], pitch forcing coefficent 
Cn_omega = 0; % [-], pitch damping coefficent 

%% Canards
n_canards = 2;
root_chord = 0.10;
spine_width = 0.015;
radial_offset = 0.02;
sweep = deg2rad(74);
spine_effectiveness = 0.5;
area_spine_eff = spine_effectiveness * spine_width * root_chord;
area_delta = 0.5 * root_chord^2 / tan(sweep);
area_canards = n_canards * (area_delta + area_spine_eff);
cp_offset = (area_delta * (root_chord / tan(sweep) / 3) - area_spine_eff * spine_width / 2) / (area_delta + area_spine_eff);
arm_canard = 0.203 / 2 + radial_offset + cp_offset;
c_canard = area_canards * arm_canard;

%% Environment
g = [-9.81; 0; 0]; % [m/s^2], gravitational acceleration in the geographic inertial frame

%% Sensors
%%% S: [1], rotation transform from sensor frame to body frame
%%% d: [m], center of sensor frame relative to body frame

S_board_acc = [0, 0,-1;
               1, 0, 0;
               0,-1, 0]; % Onboard STM IMU
S_board_mag = [0, 0,-1;
               0,-1, 0;
               1, 0, 0]; % Onboard STM IMU

S_mti = [0, 0, 1;
         1, 0, 0;
         0, 1, 0]; % Movella MTi

S_ad = [0, 0, 1;
        0,-1, 0;
        1, 0, 0]; % AD breakout board

d_board = [1.73; -0.03; -0.01]; % Onboard STM IMU
d_mti = [1.74; 0; 0]; % Movella MTi
d_ad = [1.74; -0.03; 0]; % AD breakout board

%% save and export
save("gnc/model_params.mat");