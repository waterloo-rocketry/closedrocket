function [state, cov_norm, airdata, roll_state] = navigation_codegen_entry(dt, flight_phase, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag)
    %#codegen
    % Calls the pad and flight filters.

    persistent x P b

    %% initialize at beginning
    xhat = zeros(11,1); 
    xhat(1) = 1; 
    Phat = zeros(11);

    if isempty(x)
        x = xhat; 
        P = Phat;
    end
    

    %% Pad filter iteration
    if flight_phase == false || isempty(b) % only before ignition (or if not run before)
        [xhat, bias] = pad_filter(board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag);
        x = xhat; b = bias;
    end 

    %% Flight filter iteration
    if flight_phase == true % in flight
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