function [x, P] = flight_filter(dt, x, P, bias, board_imu, mti_imu, ad_imu, board_baro, board_mag, mti_baro, mti_mag)
    % Inputs: estimates x, P; control u; measurement y; sensor bias b; timecode t
    % Input parameters: weighting Q, R; time difference to last compute T; 
    % Outputs: new estimates x, P
    %#codegen

    %% IMU Prediction + Correction steps
    %%% x = [   q(4),           w(3),         v(3),    alt(1)]
    %%% Q is a square 11 matrix, tuning for prediction E(noise)
    %%% R is a square 3 matrix, tuning for measurement E(noise) of the gyroscope
    Q = diag([[1,1,1,1]*1e-7, [10, 10, 10], [1,1,1]*1e-3, 1]) * 1e-3;
    R = diag([1, 1, 1])*1e-6;
    
    % board_imu is struct containing accel, accel_status, gyro, gyro_status, gyro_bias
    [a, w] = ekf_prefilter_imu(board_imu, mti_imu, ad_imu);
    [xhat, Phat] = ekf_update(dt, x, P, a, w, Q, R);
    x = xhat; P = Phat;

    %% Correction steps, sequential for each absolute sensor
    %%% Correction with barometer, magnetometer
    %%% R is a square matrix (size depending on amount of sensors), tuning for measurement E(noise)
 
    %%% Barometer
    if board_baro.baro_status == 1 % only correct with alive IMUs
        %%% y = [ P(1) ]
        R = 1;
        [xhat, Phat] = ekf_correct(@meas_baro, @meas_baro_jacobian, x, P, board_baro.baro_meas, bias, R);
        x = xhat; P = Phat;
    end
    if mti_baro.baro_status == 1 % only correct with alive IMUs
        %%% y = [ P(1) ]
        R = 1;
        [xhat, Phat] = ekf_correct(@meas_baro, @meas_baro_jacobian, x, P, mti_baro.baro_meas, bias, R);
        x = xhat; P = Phat;
    end

    %%% Magnetometer
    if board_mag.mag_status == 1
        %%% y = [  Mag(3) ]
        R = diag([1, 1, 1])*0.01;
        [xhat, Phat] = ekf_correct(@meas_mag, @meas_mag_jacobian, x, P, board_mag.mag_meas, bias, R);
        x = xhat; P = Phat;
    end 
    if mti_mag.mag_status == 1
        %%% y = [  Mag(3) ]
        R = diag([1, 1, 1])*0.01;
        [xhat, Phat] = ekf_correct(@meas_mag, @meas_mag_jacobian, x, P, mti_mag.mag_meas, bias, R);
        x = xhat; P = Phat;
    end 

end