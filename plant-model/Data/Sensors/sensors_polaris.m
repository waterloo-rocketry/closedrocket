%% Notes
%%% Units are m/s^2, rad/s, Gauss, Pa
% Convert from datasheet if necessary

%%% In datasheets, noise is provided as noise density in unit/sqrt(Hz), or RMS.
% Simulink White Noise needs it as height of power spectral density, so
% input the noise value as (noise density)^2, or (RMS)^2/bandwidth = (RMS)^2*samplingrate.

%% Processor values
samplingtime_nav = 0.0025; % sampling period of navigation
samplingtime_ctrl = 0.01; % sampling period of controller
samplingtime_imu = 0.0025; % sampling period imu
samplingtime_other = 0.01; % sampling period baro and mag
samplingtime_motor = 0.005;

%% Sensor group 1 (Onboard)

imu1_accel_limit = 32*9.81;
imu1_accel_bias = 12e-3*9.81; 
imu1_accel_noise = (130e-6*9.81)^2; 

imu1_gyro_limit = deg2rad(4000);
imu1_gyro_bias = deg2rad(1);
imu1_gyro_noise = (deg2rad(2.8e-3))^2;

imu1_mag_limit = 49.152; 
imu1_mag_noise = (3e-3)^2 * samplingtime_other;

imu1_baro_limit = [1000, 120e3]; % not yet used in Sim
imu1_baro_bias = 150;
imu1_baro_noise = 250;


%% Sensor group 2 values (Movella MTi)

imu2_accel_limit = 10*9.81;
imu2_accel_bias = 15e-6*9.81;
imu2_accel_noise = (60e-6*9.81)^2;

imu2_gyro_limit = deg2rad(2000);
imu2_gyro_bias = deg2rad(8/3600);
imu2_gyro_noise = (deg2rad(0.007))^2;

imu2_mag_limit = 8;
imu2_mag_noise = (1e-3)^2 * samplingtime_other;

imu2_baro_limit = [30e3, 125e3]; % not yet used in Sim
imu2_baro_bias = 8;
imu2_baro_noise = (1.2)^2 * samplingtime_other;


%% Sensor group 3 (AD breakout)

imu3_accel_limit = 16*9.81;
imu3_accel_bias = 60e-3*9.81;
imu3_accel_noise = (27e-6*9.81)^2; 

imu3_gyro_limit = deg2rad(20000);
imu3_gyro_bias = deg2rad(2.5/100);
imu3_gyro_noise = (deg2rad(0.25))^2;