function [x_init, bias] = pad_filter(board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag)
    % Computes inital state for flight filter, and bias values for the sensors
    % Outputs: initial state, sensor biases
    %#codegen

    %% Settings
    alpha = 0.0005; % low pass time constant

    %% parameters
    persistent param
    if isempty(param)
        param = coder.load("model_params.mat"); % only required parameter is launch altitude
    end 

    %% Initialization
    %%% remember filtered values from last iteration
    persistent board_accel_f board_gyro_f mti_accel_f mti_gyro_f ad_accel_f ad_gyro_f 
    persistent board_baro_f board_mag_f mti_baro_f mti_mag_f

    if isempty(board_accel_f)
        board_accel_f = board_accel.meas;
        board_gyro_f = board_gyro.meas;
        mti_accel_f = mti_accel.meas;
        mti_gyro_f = mti_gyro.meas; 
        ad_accel_f = ad_accel.meas; 
        ad_gyro_f = ad_gyro.meas;  
        board_baro_f = board_baro.meas; 
        board_mag_f = board_mag.meas; 
        mti_baro_f = mti_baro.meas; 
        mti_mag_f = mti_mag.meas; 
        bias.board_gyro = zeros(3,1); 
        bias.mti_gyro = zeros(3,1); 
        bias.ad_gyro = zeros(3,1); 
        bias.board_mag_earth = zeros(3,1); 
        bias.mti_mag_earth = zeros(3,1); 
        bias.board_baro = 0;
        bias.mti_baro = 0; 
    end

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