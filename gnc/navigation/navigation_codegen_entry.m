function [x, P, bias, sens_filt, roll_state, pdyn, w_status_nav] = navigation_codegen_entry(dt, flight_phase, x, P, bias, sens_filt, sens_in)
    %#codegen
    % Calls the pad and flight filters.

    w_status_nav = false;

    %% Pad filter iteration
    if flight_phase == false || isempty(bias) % only before ignition (or if not run before)
        [x, bias, sens_filt] = pad_filter(sens_in, sens_filt);
        w_status_nav = true;
    end

    %% Flight filter iteration
    if flight_phase == true % in flight
        [x, P] = flight_filter(dt, x, P, bias, sens_in);
        w_status_nav = false;
    end

    %% Compute air data
    airdata = airdata_atmos(x(11));
    airdata = airdata_dynamic(airdata, x(8:10));
    pdyn = airdata.dynamic_pressure;

    %% controller input vector
    phi = quaternion_to_roll(x(1:4));
    roll_state = [phi; x(5)];
end
