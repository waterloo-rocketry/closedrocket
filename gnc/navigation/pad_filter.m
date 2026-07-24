function [x_init, bias, sens_filt, status] = pad_filter(sens_in, sens_filt)
    % Computes on pad: inital state for flight filter, and bias values for the sensors
    % Outputs: initial state x, sensor biases bias
    %#codegen

    %% Settings
    alpha = 0.0005; % [s], low pass time constant


    %% parameters
    persistent param
    if isempty(param)
        param = coder.load("model_params.mat"); % only required parameter is launch altitude
    end

    %% lowpass filter
    % filtered = filtered + alpha*(measured-filtered);

    sens_filt.board_accel = pad_lowpass(sens_filt.board_accel, sens_in.board_accel, alpha);
    sens_filt.board_gyro = pad_lowpass(sens_filt.board_gyro, sens_in.board_gyro, alpha);
    sens_filt.mti_accel = pad_lowpass(sens_filt.mti_accel, sens_in.mti_accel, alpha);
    sens_filt.mti_gyro = pad_lowpass(sens_filt.mti_gyro, sens_in.mti_gyro, alpha);
    sens_filt.ad_accel = pad_lowpass(sens_filt.ad_accel, sens_in.ad_accel, alpha);
    sens_filt.ad_gyro = pad_lowpass(sens_filt.ad_gyro, sens_in.ad_gyro, alpha);
    sens_filt.board_baro = pad_lowpass(sens_filt.board_baro, sens_in.board_baro, alpha);
    sens_filt.board_mag = vec_norm(pad_lowpass(sens_filt.board_mag, sens_in.board_mag, alpha));
    sens_filt.mti_baro = pad_lowpass(sens_filt.mti_baro, sens_in.mti_baro, alpha);
    sens_filt.mti_mag = vec_norm(pad_lowpass(sens_filt.mti_mag, sens_in.mti_mag, alpha));

    %% Initial state determination
    %%% Orientation
    a = sens_filt.board_accel + sens_filt.mti_accel + sens_filt.ad_accel;
    [q, status] = pad_inclinometer(a); % a gets normed inside function

    %%% launch altitude
    alt = param.altitude_initial;

    %%% set constant initials
    w = [0; 0; 0]; % stationary on rail
    v = [0; 0; 0]; % stationary on rail

    %%% conconct state vector
    x_init = [q; w; v; alt];

    %% Bias determination

    %%% gyroscope
    bias.board_gyro = sens_filt.board_gyro;
    bias.mti_gyro = sens_filt.mti_gyro;
    bias.ad_gyro = sens_filt.ad_gyro;

    %%% earth magnetic field
    ST = transpose(quaternion_rotmatrix(q)); % launch attitude
    bias.board_mag_earth = ST * sens_filt.board_mag;
    bias.mti_mag_earth = ST * sens_filt.mti_mag;

    %%% barometer
    pressure = airdata_atmos(param.altitude_initial).pressure; % pressure at launch elevation
    bias.board_baro = sens_filt.board_baro - pressure;
    bias.mti_baro = sens_filt.mti_baro - pressure;

end
