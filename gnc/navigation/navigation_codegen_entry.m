function [x, P, bias, sens_filt, airdata, roll_state] = navigation_codegen_entry(dt, flight_phase, x, P, bias, sens_filt, sens_input)
    %#codegen
    % Calls the pad and flight filters.

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
        [x, bias, sens_filt] = pad_filter(board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag, sens_filt);
    end

    %% Flight filter iteration
    if flight_phase == true % in flight
        [x, P] = flight_filter(dt, x, P, bias, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag);
    end

    %% Compute air data
    airdata = airdata_atmos(x(11));
    airdata = airdata_dynamic(airdata, x(8:10));

    %% controller input vector
    phi = quaternion_to_roll(x(1:4));
    roll_state = [phi; x(5)];

end
