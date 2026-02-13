function [xhat, Phat_norm, airdata, xR] = navigation_module(timestamp, IMU_1, IMU_2, sensor_select)
    % Top-level estimator module. Calls EKF algorithm.
    % Inputs: concocted measurement and output vectors with multiple sensors. Not yet fully supported, work in progress
    % IMU = struct of IMUi = [accel; omega; mag; baro] 
    %#codegen
    
    persistent t x P b flight_phase; % remembers t, x, P from last iteration
    
    %% settings
    idle_time = 9; % wait time to handover

    %% initialize at beginning
    xhat = zeros(11,1); xhat(1) = 1; Phat = zeros(11); bias_1 = zeros(10, 1); bias_2 = zeros(10, 1);
    if isempty(x)
        x = xhat; P = Phat; b.bias_1 = bias_1; b.bias_2 = bias_2;
        flight_phase = 1;
        if timestamp >= 0.005
                t = timestamp-0.005;
        else 
                t = 0;
        end
    end
    
    %% timecode
    dt = timestamp - t; % time step size for integrators
    t = timestamp;
    
    %% flight phase
    if t >= idle_time % mock for flight phase
            flight_phase = 0; % 1 is pad, 0 is flight
    end

    %% Pad filter iteration
    if flight_phase ~= 0 % only before ignition
        [xhat, bias_1, bias_2] = pad_filter(IMU_1, IMU_2, sensor_select);
        x = xhat; b.bias_1 = bias_1; b.bias_2 = bias_2;
    end 

    %% Flight filter iteration
    if flight_phase == 0 % in flight
        [xhat, Phat] = flight_filter(dt, x, P, b, IMU_i);
        x = xhat; P = Phat;
    end

    %% Compute air data
    airdata = airdata_dynamic(alt, v);

    %% Compute variance norm for EKF quality
    Phat_norm = norm(P); % Compute the norm of the covariance matrix

    %% Roll state for controller
    phi = quaternion_to_roll(x(1:4));
    xR = [phi; x(5)];
end