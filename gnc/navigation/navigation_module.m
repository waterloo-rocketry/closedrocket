function [xhat, Phat_norm, airdata, xR] = navigation_module(timestamp, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag)
    % Top-level estimator module. Calls the pad and flight filters.
    %#codegen    
    persistent t x P b flight_phase; % remembers t, x, P from last iteration
    
    %% settings
    idle_time = 9; % wait time to handover

    %% initialize at beginning
    xhat = zeros(11,1); xhat(1) = 1; Phat = zeros(11); 
    if isempty(x)
        x = xhat; P = Phat;
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
    if flight_phase ~= 0 || isempty(b) % only before ignition (or if not run before)
        [xhat, bias] = pad_filter(board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag);
        x = xhat; b = bias;
    end 

    %% Flight filter iteration
    if flight_phase == 0 % in flight
        [xhat, Phat] = flight_filter(dt, x, P, b, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag);
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