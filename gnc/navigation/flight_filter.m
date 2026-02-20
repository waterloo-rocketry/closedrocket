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
    R = diag([[1, 1, 1])*1e-6;
    
    A1 = IMU_1(1:3,1);
    W1 = IMU_1(4:6,1);
    A2 = IMU_2(1:3,1);
    W2 = IMU_2(4:6,1);

    [a, w] = ekf_prefilter_imu(dt, bias, A1, W1, A2, W2, A3, W3);
    [xhat, Phat] = ekf_update(dt, x, P, a, w, Q, R);
    x = xhat; P = Phat;

    %% Correction steps, sequential for each absolute sensor
    %%% Correction with barometer, magnetometer
    %%% R is a square matrix (size depending on amount of sensors), tuning for measurement E(noise)
 
    % Barometer
    if sensor_select(1) == 1 % only correct with alive IMUs
        %%% y = [ P(1) ]
        R = 1;

        [xhat, Phat] = ekf_correct(@meas_baro, @meas_baro_jacobian, x, P, IMU_1(10), bias, R);
        x = xhat; P = Phat;
    end


    % Magnetometer
    if sensor_select(2) == 1
        %%% y = [  Mag(3) ]
        R = diag([1, 1, 1])*0.01;

        [xhat, Phat] = ekf_correct(@meas_mag, @meas_mag_jacobian, x, P, IMU_2(7:9), bias, R);
        x = xhat; P = Phat;
    end 

end