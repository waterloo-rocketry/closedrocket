function [x, P, bias, sens_filt, cov_norm, roll_state, pdyn] = navigation_codegen_entry(dt, flight_phase, x, P, bias, sens_filt, sens_in)
    %#codegen
    % Calls the pad and flight filters.

    %% Pad filter iteration
    if flight_phase == false || isempty(bias) % only before ignition (or if not run before)
        [x, bias, sens_filt] = pad_filter(sens_in, sens_filt);
    end

    %% Flight filter iteration
    if flight_phase == true % in flight
        [x, P] = flight_filter(dt, x, P, bias, sens_in);
    end

    %% Compute air data
    airdata = airdata_atmos(x(11));
    airdata = airdata_dynamic(airdata, x(8:10));
    pdyn = airdata.dynamic_pressure;

    %% Compute variance norm
    %%% for evaluating EKF quality
    cov_norm = norm(P, inf); % scalar, infinity norm of the covariance matrix

    %% controller input vector
    phi = quaternion_to_roll(x(1:4));
    roll_state = [phi; x(5)];

end
