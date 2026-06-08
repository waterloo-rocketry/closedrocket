%% Data preparation for simulation
or_data = fillmissing(or_data, 'previous');

%% Constant mass and inertia
% Wet
I_xx_0 = or_data.RotationalMomentOfInertia_kg_m__(1);
I_yy_0 = or_data.LongitudinalMomentOfInertia_kg_m__(1);
I_zz_0 = I_yy_0; % rocket is axially symmetric
I_0 = [I_xx_0 0 0; 0 I_yy_0 0; 0 0 I_zz_0]; %kg*m^2
total_mass_0 = or_data.Mass_g_(1) / 1000; %g -> kg

% Dry
I_xx_d = or_data.RotationalMomentOfInertia_kg_m__(find(~isnan(or_data.RotationalMomentOfInertia_kg_m__), 1, 'last')); %moments of inertia become NaN at the end for some reason?
I_yy_d = or_data.LongitudinalMomentOfInertia_kg_m__(find(~isnan(or_data.LongitudinalMomentOfInertia_kg_m__), 1, 'last'));
I_zz_d = I_yy_d; % rocket is axially symmetric
I_d = [I_xx_d 0 0; 0 I_yy_d 0; 0 0 I_zz_d]; %kg*m^2
total_mass_d = or_data.Mass_g_(find(~isnan(or_data.Mass_g_), 1, 'last')) / 1000; % g -> kg

%% Input Time Series
sim_time = or_data.x_Time_s_; % s

%%% time series mass and inertia
mass = or_data.Mass_g_ / 1000; % g -> kg
OR_cg = -or_data.CGLocation_cm_ / 100; % centre of gravity over time (replace NaN with Last val)
I_xx = or_data.RotationalMomentOfInertia_kg_m__;
I_yy = or_data.LongitudinalMomentOfInertia_kg_m__;
I_zz = I_yy; % rocket is axially symmetric

%%% thrust
F_thrust = or_data.Thrust_N_; % N assume thrust is perfectly aligned
F_thrust = circshift(F_thrust, 5);
F_thrust(end-5:end) = 0;
F_thrust(2:5) = NaN;
F_thrust = fillmissing(F_thrust, 'linear');

%%% reference trajectory 
alt = or_data.Altitude_m_;
vel_vert = or_data.VerticalVelocity_m_s_;
drag = or_data.DragForce_N_;
crossrange = or_data.LateralDistance_m_;
v_crossrange = or_data.LateralVelocity_m_s_;
angle_vert_deg = 90 - or_data.VerticalOrientation_zenith____;

%%% import Mach/Cd
cd_time = or_data.DragCoefficient___;
mach_time = or_data.MachNumber___;

%% Cd vs Mach table
% Combine the two vectors into a 2-column matrix
cd_data_raw = [mach_time(:), cd_time(:), vel_vert(:)];
% Remove rows where vertical velocity is negative
cd_neg_vel_idx = cd_data_raw(:,3) <= 0; % logical index for rows with vel_vert < 0
cd_data_raw(cd_neg_vel_idx, :) = []; % remove those rows from cd_data

% Remove rows with duplicate Mach values (keep the first occurrence)
[~, cd_unique_idx] = unique(cd_data_raw(:, 1), 'first');  % Find unique Mach numbers
cd_data_unique = cd_data_raw;  % Keep only rows with unique Mach numbers

% Remove rows with any NaN values
cd_data_unique = cd_data_unique(~any(isnan(cd_data_unique), 2), :);  % Remove rows with NaN

% Sort the rows by the Mach number (first column)
cd_data_sorted = sortrows(cd_data_unique, 1);  % Sort by Mach number (first column)

% Get the Cd value corresponding to Mach 0.3
mach_01_idx = find(cd_data_sorted(:, 1) >= 0.1, 1, 'first');  % Find the index of Mach closest to or equal to 0.3
cd_at_01 = cd_data_sorted(mach_01_idx, 2);  % Cd value at Mach 0.3

% Replace Cd values for Mach < 0.3 with Cd value at Mach 0.3
cd_data_sorted(cd_data_sorted(:, 1) < 0.3, 2) = cd_at_05;

% smooth (necessary because using OR export is dumb)
cd_data_smooth = cd_data_sorted;
cd_data_smooth(:,2) = movmean(cd_data_smooth(:,2), 40);

% Final lookup table
cd_lookup_table = cd_data_smooth(:,1:2);
cd_input = cd_lookup_table(:, 1); % Mach #
cd_data = cd_lookup_table(:, 2); % Cd


%% Aero scripts
% Nose
[nosecone_area_planform, nosecone_area_bow, nosecone_area_aft, nosecone_volume, nosecone_pos_x_cp] = aerobody(nosecone_length, 0, 2 * nosecone_radius, 0);


% Body
[body_area_planform, body_area_bow, body_area_aft, body_volume, body_pos_x_cp] = aerobody(body_length, rocket_diameter, rocket_diameter, - nosecone_length);

% Fins 
% [fin_pos_x_cp, fin_Cnfdelta, fin_CndNi, fin_CNa, fin_aspectratio, fin_area, fin_midchord_angle, fin_dist_chord_mean, fin_pos_r_chord_mean, fin_leading_edge] = fins(fin_chord_root, fin_chord_tip, fin_height, fin_sweep, fin_pos_x_roottip, fin_number, rocket_area_frontal, rocket_diameter);
[fin_pos_x_cp, fin_pos_r_mean, fin_area, fin_aspectratio, fin_midchord_angle, fin_factor, fin_pos_x_cp_mach2] = aerosurface(fin_chord_root, fin_chord_tip, fin_height, fin_sweep_angle, fin_pos_x_roottip, fin_number, rocket_diameter);


% Tail
[tail_area_planform, tail_area_bow, tail_area_aft, tail_volume, tail_pos_x_cp] = aerobody(tail_length, 2 * tail_radius_outer, 2 * tail_radius_smallest, - nosecone_length - body_length);

% Canards 
% [canard_pos_x, canard_Cnalfat, canard_Cnfdelta, canard_CndNi, canard_aspectratio, canard_area, canard_midchord_angle, canard_dist_chord_mean, canard_pos_r_mean, canard_leading_edge] = canards(canard_chord_root, canard_chord_tip, canard_height, canard_pos_x_roottip, canard_number, rocket_area_frontal, rocket_diameter);
[canard_pos_x, canard_pos_r_mean, canard_area, canard_aspectratio, canard_midchord_angle, canard_factor, canard_pos_x_cp_mach2] = aerosurface(canard_chord_root, canard_chord_tip, canard_height, canard_sweep_angle, canard_pos_x_roottip, canard_number, rocket_diameter);