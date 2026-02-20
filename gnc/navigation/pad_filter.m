function [x_init, bias] = pad_filter(board_imu, mti_imu, ad_imu, board_baro, board_mag, mti_baro, mti_mag)
    % Computes inital state and covariance estimate for EKF, and bias values for the IMU
    % Uses all available sensors: Gyroscope W, Magnetometer M, Accelerometer A, Barometer P
    % Outputs: initial state, sensor bias matrix
    %#codegen

    %%% remember filtered values from last iteration
    persistent board_accel_f board_gyro_f mti_accel_f mti_gyro_f ad_accel_f ad_gyro_f 
    persistent board_baro_f board_mag_f mti_baro_f mti_mag_f

    %% parameters
    persistent param
    if isempty(param)
        param = coder.load("model/model_params.mat");
    end 

    %% lowpass filter (and initialization)
    alpha = 0.0005; % low pass time constant
    % filtered = filtered + alpha*(measured-filtered);

    if board_imu.accel_status == 1, board_accel_f = lowpass(board_accel_f, board_imu.accel_meas, alpha); end
    if board_imu.gyro_status == 1, board_gyro_f = lowpass(board_gyro_f, board_imu.gyro_meas, alpha); end
    if mti_imu.accel_status == 1, mti_accel_f = lowpass(mti_accel_f, mti_imu.accel_meas, alpha); end
    if mti_imu.gyro_status == 1, mti_gyro_f = lowpass(mti_gyro_f, mti_imu.gyro_meas, alpha); end
    if ad_imu.accel_status == 1, ad_accel_f = lowpass(ad_accel_f, ad_imu.accel_meas, alpha); end
    if ad_imu.gyro_status == 1, ad_gyro_f = lowpass(ad_gyro_f, ad_imu.gyro_meas, alpha); end
    if board_baro.baro_status == 1, board_baro_f = lowpass(board_baro_f, board_baro.baro_meas, alpha); end
    if board_mag.mag_status == 1, board_mag_f = lowpass(board_mag_f, board_mag.mag_meas, alpha); end
    if mti_baro.baro_status == 1, mti_baro_f = lowpass(mti_baro_f, mti_baro.baro_meas, alpha); end
    if mti_mag.mag_status == 1, mti_mag_f = lowpass(mti_mag_f, mti_mag.mag_meas, alpha); end


    %% State determination
    
    %%% average specific force of selected sensors
    a = zeros(3,1); % acceleration 
    if board_imu.accel_status == 1 % only add alive IMUs to average
        a = a + board_accel_f;
    end
    if mti_imu.accel_status == 1
        a = a + mti_accel_f;
    end
    if ad_imu.accel_status == 1
        a = a + ad_accel_f;
    end

    %%% gravity vector in body-fixed frame
    A = a / (norm(a) + 1e-6); % unit vector of gravity direction
    
    %%% determine initial orientation quaternion
    qw = sqrt( 0.5 + 0.5*A(1) );
    qx = 0;
    if qw == 0 % exact upside down case
        qy = 1; % either qy = 1 or qz = 1, this is arbitrary 
        qz = 0;
    else 
        qy = 0.5 * A(3) / qw;
        qz = -0.5 * A(2) / qw;
    end
    q = [qw; qx; qy; qz];
    q = q / norm(q);

    %%% launch altitude
    alt = param.elevation;
    
    %%% set constant initials
    w = [0; 0; 0]; % stationary on rail
    v = [0; 0; 0]; % stationary on rail

    %%% conconct state vector
    x_init = [q; w; v; alt];

    %% Bias determination
        
    %%% gyroscope
    if board_imu.gyro_status == 1
        bias.board_gyro = board_gyro_f;
    end
    if mti_imu.gyro_status == 1
        bias.mti_gyro = mti_gyro_f;
    end
    if ad_imu.gyro_status == 1
        bias.ad_gyro = ad_gyro_f;
    end
    
    %%% earth magnetic field
    ST = transpose(quaternion_rotmatrix(q)); % launch attitude
    if board_imu.gyro_status == 1
        bias.board_mag_earth = ST * board_mag_f;
    end
    if mti_imu.gyro_status == 1
        bias.mti_mag_earth = ST * mti_mag_f;
    end

    %%% barometer
    pressure = model_airdata(param.elevation).pressure; % what the pressure should be at launch elevation
    if board_baro.baro_status == 1
        bias.board_baro = board_baro_f - pressure;
    end
    if mti_baro.baro_status == 1
        bias.mti_baro = mti_baro_f - pressure;
    end

end

%% Lowpass filter function
function filtered = lowpass(filtered, measured, alpha)
    if isempty(filtered) % initialize
        filtered = measured;  
    else 
        filtered = alpha * measured + (1 - alpha) * filtered; 
    end
end