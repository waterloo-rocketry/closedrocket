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
S1 = [0, 0, 1;
     1, 0, 0;
     0, 1, 0]; % IMU 1, rotation transform from sensor frame to body frame
S2 = [0, 0, -1;
     -1, 0, 0;
      0, 1, 0]; % IMU 2, rotation transform from sensor frame to body frame

d1 = [1.2; 0.074; -0.027]; % center of sensor frame
d2 = [1.2; 0.065; 0.047]; % center of sensor frame

B1 = eye(3); % Soft iron compensation
B2 = eye(3); % Soft iron compensation
b1 = [0;0;0]; % Hard iron compensation
b2 = [0;0;0]; % Hard iron compensation

%% save and export
save("design/model/model_params.mat");