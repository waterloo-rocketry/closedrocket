function [x_ret, P_ret, b_ret, sf_ret] = navigation_codegen_entry(dt, flight_phase, x, P, b, sf, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag)
    %#codegen
    % Calls the pad and flight filters.
    
    x_ret = x;
    P_ret = P;
    b_ret = b;
    sf_ret = sf;

    %% Pad filter iteration
    if flight_phase == false || isempty(b) % only before ignition (or if not run before)
        [xhat, bias, sf_ret] = pad_filter(board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag, sf);
        x_ret = xhat; b_ret = bias;
    end 

    %% Flight filter iteration
    if flight_phase == true % in flight
        [x_ret, P_ret] = flight_filter(dt, x, P, b, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag);
    end

end