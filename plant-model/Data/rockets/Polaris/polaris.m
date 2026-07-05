%% OR Simulation Output Data
or_data = readtable('plant-model/data/rockets/Polaris/polaris_cycle_2.csv');

%% Initial values
location = [420; 43.47; -80.54]; % launch location on earth. Altitude, Latitude, Longitude [m, deg, deg]
rail_angle_pitch = deg2rad(-6); % Rail pitch angle. Negative is pitched downrange [rad]
rail_angle_yaw = deg2rad(0); % Rail yaw angle. Negative is yawed downrange [rad]
rail_angle_roll = deg2rad(0); % Rocket clocking angle [rad]
rail_length = 8.23; % [m]

%% Sensor mounting
%%% _d : Mounting location relative to body frame [m]
%%% _S : Mounting orientation relative to nosetip [rotation matrix]
%%%      From body to sensor
%%% Sensor group 1 (Onboard)
sensor_1_d = [-1.51; -0.03; -0.01]; 
sensor_1_S_acc = [0, 1, 0;
                  0, 0,-1;
                 -1, 0, 0]; 
sensor_1_S_mag = [0, 0, 1;
                  0,-1, 0;
                 -1, 0, 0]; 
%%% Sensor group 2 (Movella MTi)
sensor_2_d = [-1.52; 0.0; 0.0];
sensor_2_S = [0, 1, 0;
              0, 0, 1;
              1, 0, 0]; 

%%% Sensor group 3 (AD breakout)
sensor_3_d = [-1.52; -0.03; 0.0]; 
sensor_3_S = [0, 0, 1;
              0,-1, 0;
              1, 0, 0];
%% Actuator parameters
act_freq = 167; % natural frequency, approx 1/timeconstant [1/s]
act_deadtime = 0.003; % delay in servo internal control loop [s]
act_damping = 0.9; % damping ratio
act_backlash = 0.3; % play [deg]
act_anglelimit = 22; % max deflection [deg]
act_ratelimit = 600; % max rate [deg/s]
act_gear_ratio = - 0.5; % speed ratio of gearing between motor and canard

%% Misc Rocket parameters
engine_thrust_factor = 1; % perfomance gain
drag_factor = 1; % drag gain
canard_roll_reversal_factor = 1; % coefficient gain

%% Aerodynamics Reference Geometry
%Reference parameters   
rocket_diameter = 0.203; % reference length [m]
rocket_area_frontal = pi * rocket_diameter^2 / 4; % reference area [m^2]
rocket_length = 5.58; % rocket length [m]

%Parachutes
chute_pos_x = -1.1; % chute attachment [m]
chute_drogue_drag = 0.55 * 0.82; % Cd * A [m^2]
chute_drogue_time = 2; % time after apogee [s]
chute_main_drag = 1.23 * 8.36; % Cd * A [m^2]
chute_main_time = 0; % time after threshold altitude [s]
chute_main_threshold = 2000; % threshold altitude [m]

%Nosecone parameters
nosecone_length = 1.02; % nosecone length [m]
nosecone_radius = rocket_diameter / 2; % nosecone radius [m]

%Tail parameters
tail_radius_outer = 0.203 / 2; % tail radius [m]
tail_length = 0.0254; % tail length [m]
tail_radius_smallest = 0.19 / 2; % smallest tail radius(?) [m]
tail_pos_x_roottip = -rocket_length + tail_length; % tail position measured from nosecone

%Body parameters
body_length = rocket_length - nosecone_length - tail_length; % fuselage length only [m]
body_surface_roughness = 20 / 10^6; % RMC(?) roughness 20 um smooth paint

%Fin parameters
fin_chord_root = 0.483; %[m] root chord?
fin_chord_tip = 0.114; %[m] tip chord?
fin_height = 0.229; %[m] height?
fin_sweep = 0.432; % [m]
fin_sweep_angle = deg2rad(62.1); % angle from radial normal [rad]
fin_pos_x_roottip = - ( rocket_length - tail_length - 0.483 - 0.0127); % position of fins measured from nosecone [m]
fin_number = 4; % Number of fins
fin_cant_angle_rad = deg2rad(0); % fin cant angle [rad]

% Canards parameters 
canard_number = 2;
canard_chord_root = 0.10; % root chord
canard_chord_tip = 0; % tip chord
canard_sweep_angle = deg2rad(74); % angle from radial normal [rad]
canard_height = canard_chord_root / tan(canard_sweep_angle); % root to tip length
canard_spine_width = 0.015;
canard_spine_effectiveness = 0.5;
canard_radial_offset = 0.02;
canard_delta_max = deg2rad(12); % Canard maximum angle of attack
canard_pos_x_roottip = - (nosecone_length + 0.508 + 0.102 / 2 - 0.03); % approximate position of the most forward tip of the canards
canard_cant_zero = deg2rad(0.2); % zero roll not perfect