function [x_init, bias, sens_filt] = pad_filter(board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag, sens_filt)
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

    sens_filt.board_accel = pad_lowpass(sens_filt.board_accel, board_accel, alpha);
    sens_filt.board_gyro = pad_lowpass(sens_filt.board_gyro, board_gyro, alpha);
    sens_filt.mti_accel = pad_lowpass(sens_filt.mti_accel, mti_accel, alpha);
    sens_filt.mti_gyro = pad_lowpass(sens_filt.mti_gyro, mti_gyro, alpha);
    sens_filt.ad_accel = pad_lowpass(sens_filt.ad_accel, ad_accel, alpha);
    sens_filt.ad_gyro = pad_lowpass(sens_filt.ad_gyro, ad_gyro, alpha);
    sens_filt.board_baro = pad_lowpass(sens_filt.board_baro, board_baro, alpha);
    sens_filt.board_mag = pad_lowpass(sens_filt.board_mag, board_mag, alpha);
    sens_filt.mti_baro = pad_lowpass(sens_filt.mti_baro, mti_baro, alpha);
    sens_filt.mti_mag = pad_lowpass(sens_filt.mti_mag, mti_mag, alpha);

    %% Initial state determination
    %%% Orientation
    a = zeros(3,1); % acceleration
    if board_accel.status == 1 % only add alive IMUs to average
        a = a + sens_filt.board_accel;
    end
    if mti_accel.status == 1
        a = a + sens_filt.mti_accel;
    end
    if ad_accel.status == 1
        a = a + sens_filt.ad_accel;
    end
    q = pad_inclinometer(a); % a gets normed inside function

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
    pressure = airdata_atmos(param.elevation).pressure; % what the pressure should be at launch elevation
    bias.board_baro = sens_filt.board_baro - pressure;
    bias.mti_baro = sens_filt.mti_baro - pressure;

end
