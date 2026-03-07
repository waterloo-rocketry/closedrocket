function [xhat, Phat_norm, airdata, xR] = navigation_module(timestamp, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag)
    % Top-level estimator module. Calls the pad and flight filters.
    %#codegen    
    persistent t x P b flight_phase k; % remembers t, x, P from last iteration
    
    %% settings
    idle_time = 9; % wait time to handover

    %% initialize at beginning
    xhat = zeros(11,1); xhat(1) = 1; Phat = zeros(11); 
    if isempty(x)
        x = xhat; P = Phat;
        flight_phase = 1;
        k = 1;
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

    %% Mock slower sampling rate of baro, mag
    sampling_imu = 0.002; % sampling period of imu [s]
    sampling_other = 0.02; % sampling period of baro, mag [s]
    if k ==  sampling_other/sampling_imu
        k = 1; % reset counter for the next cycle
    else
        k = k + 1; % increment counter for the next cycle
        board_baro.status = false; % sensors do not have new data
        board_mag.status = false;
        mti_baro.status = false;
        mti_mag.status = false;
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
    airdata = airdata_dynamic(x(11), x(8:10));

    %% Compute variance norm for EKF quality
    Phat_norm = norm(P); % Compute the norm of the covariance matrix

    %% Roll state for controller
    phi = quaternion_to_roll(x(1:4));
    xR = [phi; x(5)];
end