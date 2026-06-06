function [x_ret, P_ret, bias_ret, sens_filt_ret, cov_norm, airdata, roll_state] = navigation_codegen_entry(dt, flight_phase, x, P, bias, sens_filt, sens_input)
    %#codegen
    % Calls the pad and flight filters.
    
    x_ret = x;
    P_ret = P;
    bias_ret = bias;
    sens_filt_ret = sens_filt;
    board_accel = sens_input.board_accel;
    board_gyro = sens_input.board_gyro;
    mti_accel = sens_input.mti_accel;
    mti_gyro = sens_input.mti_gyro;
    ad_accel = sens_input.ad_accel;
    ad_gyro = sens_input.ad_gyro;
    board_baro = sens_input.board_baro;
    board_mag = sens_input.board_mag;
    mti_baro = sens_input.mti_baro;
    mti_mag = sens_input.mti_mag;

    %% Pad filter iteration
    if flight_phase == false || isempty(bias) % only before ignition (or if not run before)
        [xhat, bias, sens_filt_ret] = pad_filter(board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag, sens_filt);
        x_ret = xhat; bias_ret = bias;
    end 

    %% Flight filter iteration
    if flight_phase == true % in flight
        [x_ret, P_ret] = flight_filter(dt, x, P, bias, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag);
    end

    %% return the other values
     %% Compute variance norm 
    %%% for evaluating EKF quality
    cov_norm = norm(P_ret); % scalar, 2-norm of the covariance matrix

    %% Compute air data
    airdata = airdata_atmos(x_ret(11));
    airdata = airdata_dynamic(airdata, x_ret(8:10));

    %% controller input vector
    phi = quaternion_to_roll(x_ret(1:4));
    roll_state = [phi; x_ret(5)];

end