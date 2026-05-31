function [x_init, bias, sf_ret] = pad_filter(board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag, sf)
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

    board_accel_f = sf.board_accel_f;
    board_gyro_f = sf.board_gyro_f;
    mti_accel_f = sf.mti_accel_f;
    mti_gyro_f = sf.mti_gyro_f;
    ad_accel_f = sf.ad_accel_f;
    ad_gyro_f = sf.ad_gyro_f;
    board_baro_f = sf.board_baro_f;
    board_mag_f = sf.board_mag_f;
    mti_baro_f = sf.mti_baro_f;
    mti_mag_f = sf.mti_mag_f;


    %% lowpass filter
    % filtered = filtered + alpha*(measured-filtered);
        
    board_accel_f = pad_lowpass(board_accel_f, board_accel, alpha);
    board_gyro_f = pad_lowpass(board_gyro_f, board_gyro, alpha); 
    mti_accel_f = pad_lowpass(mti_accel_f, mti_accel, alpha);
    mti_gyro_f = pad_lowpass(mti_gyro_f, mti_gyro, alpha); 
    ad_accel_f = pad_lowpass(ad_accel_f, ad_accel, alpha);
    ad_gyro_f = pad_lowpass(ad_gyro_f, ad_gyro, alpha); 
    board_baro_f = pad_lowpass(board_baro_f, board_baro, alpha); 
    board_mag_f = pad_lowpass(board_mag_f, board_mag, alpha); 
    mti_baro_f = pad_lowpass(mti_baro_f, mti_baro, alpha); 
    mti_mag_f = pad_lowpass(mti_mag_f, mti_mag, alpha);

    sf_ret.board_accel_f = board_accel_f;
    sf_ret.board_gyro_f = board_gyro_f;
    sf_ret.mti_accel_f = mti_accel_f;
    sf_ret.mti_gyro_f = mti_gyro_f;
    sf_ret.ad_accel_f = ad_accel_f;
    sf_ret.ad_gyro_f = ad_gyro_f;
    sf_ret.board_baro_f = board_baro_f;
    sf_ret.board_mag_f = board_mag_f;
    sf_ret.mti_baro_f = mti_baro_f;
    sf_ret.mti_mag_f = mti_mag_f;


    %% Initial state determination    
    %%% Orientation
    a = zeros(3,1); % acceleration 
    if board_accel.status == 1 % only add alive IMUs to average
        a = a + board_accel_f;
    end
    if mti_accel.status == 1
        a = a + mti_accel_f;
    end
    if ad_accel.status == 1
        a = a + ad_accel_f;
    end
    q = pad_inclinometer(a); % a gets normed inside function

    %%% launch altitude
    alt = param.elevation;
    
    %%% set constant initials
    w = [0; 0; 0]; % stationary on rail
    v = [0; 0; 0]; % stationary on rail

    %%% conconct state vector
    x_init = [q; w; v; alt];

    %% Bias determination
        
    %%% gyroscope
    bias.board_gyro = board_gyro_f;
    bias.mti_gyro = mti_gyro_f;
    bias.ad_gyro = ad_gyro_f;
    
    %%% earth magnetic field
    ST = transpose(quaternion_rotmatrix(q)); % launch attitude
    bias.board_mag_earth = ST * board_mag_f;
    bias.mti_mag_earth = ST * mti_mag_f;
    
    %%% barometer
    pressure = airdata_atmos(param.elevation).pressure; % what the pressure should be at launch elevation
    bias.board_baro = board_baro_f - pressure;
    bias.mti_baro = mti_baro_f - pressure;

end