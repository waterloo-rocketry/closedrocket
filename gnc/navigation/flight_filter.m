function [x, P] = flight_filter(dt, x, P, bias, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag)
    % Computes state in flight
    % Input variables: time step dt, old state x, old covariance P, 
    % Input parameters: sensor biases bias;
    % Input measurements: sensorgroup_sensortype
    % Outputs: new state x, new covariance P
    %#codegen

    %% IMU Prediction + Correction steps
    %%% x = [   q(4),           w(3),         v(3),    alt(1)]
    %%% Q is a square 11 matrix, tuning for expected dynamics noise magnitude E(noise)
    %%% R is a square 3 matrix, tuning for expected measurement noise magnitude E(noise) of the gyroscope
    Q = diag([[1,1,1,1]*1e-10, [1, 1, 1]*1e-2, [1,1,1]*1e-4, 1e-3]);
    R = diag([1, 1, 1])*1e-9;

    [a, w] = ekf_prefilter_imu(bias, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro);
    [xhat, Phat] = ekf_update(dt, x, P, a, w, Q, R);
    x = xhat; P = Phat;

    %% Correction steps, sequential for each additional sensor
    %%% Correction with barometer, magnetometer
    %%% R is a square matrix (size length of sensor vector), tuning for expected measurement noise magnitude E(noise)
 
    %%% Barometer
    if board_baro.status == true % only correct with alive IMUs
        %%% y = [ P(1) ]
        R = 100;
        [xhat, Phat] = ekf_correct(@meas_baro, @meas_baro_jacobian, x, P, board_baro.meas, bias.board_baro, R);
        x = xhat; P = Phat;
    end
    if mti_baro.status == true % only correct with alive IMUs
        %%% y = [ P(1) ]
        R = 100;
        [xhat, Phat] = ekf_correct(@meas_baro, @meas_baro_jacobian, x, P, mti_baro.meas, bias.mti_baro, R);
        x = xhat; P = Phat;
    end

    %%% Magnetometer
    if board_mag.status == true
        %%% y = [  Mag(3) ]
        R = diag([1, 1, 1]);
        [xhat, Phat] = ekf_correct(@meas_mag, @meas_mag_jacobian, x, P, board_mag.meas, bias.board_mag_earth, R);
        x = xhat; P = Phat;
    end 
    if mti_mag.status == true
        %%% y = [  Mag(3) ]
        R = diag([1, 1, 1]);
        [xhat, Phat] = ekf_correct(@meas_mag, @meas_mag_jacobian, x, P, mti_mag.meas, bias.mti_mag_earth, R);
        x = xhat; P = Phat;
    end 

end