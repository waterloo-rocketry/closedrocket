function [airdata] = airdata_dynamic(airdata, v)
    % appends airadata struct with dynamic air data from velocity vector
    % dynamic air data: dynamic pressure, mach number

    airspeed = norm(v);

    % return values
    airdata.mach = airspeed / airdata.sonic_speed; % [-], Mach number
    airdata.dynamic_pressure = 0.5 * airdata.density * airspeed^2; % [Pa]
end