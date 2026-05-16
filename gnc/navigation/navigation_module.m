function [state, cov_norm, airdata, roll_state] = navigation_module(timestamp, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag)
    % Top-level navigation module. Calls the pad and flight filters. 
    % Mocks firmware higher-level stuff
    %#codegen    
    persistent t x P b flight_phase k; % remembers t, x, P from last iteration
    
    %% config settings
    idle_time = 10; % [s], wait time to handover from pad to flight
    sampling_imu = 0.002; % [s], sampling period of imu
    sampling_other = 0.02; % [s], sampling period of baro, mag

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
    dt = timestamp - t; % [s], time step size for dynamics integration
    t = timestamp;
    
    %% flight phase
    if t >= idle_time % mock for flight phase
            flight_phase = 0; % 1 is pad, 0 is flight
    end

    %% Mock slower sampling rate of baro, mag
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

    %% Pack state as struct
    %%% use union in C or smth
    state.q = x(1:4); % [1], attitude quaternion
    state.w = x(5:7); % [rad/s], angular rate 
    state.v = x(8:10); % [m/s], velocity
    state.alt = x(11); % [m], altitude
    state.x = x; % also full state as vector is needed in simulink

    %% Compute variance norm 
    %%% for evaluating EKF quality
    cov_norm = norm(P); % scalar, 2-norm of the covariance matrix

    %% Compute air data
    airdata = airdata_atmos(x(11));
    airdata = airdata_dynamic(airdata, x(8:10));

    %% controller input vector
    phi = quaternion_to_roll(x(1:4));
    roll_state = [phi; x(5)];
end