function [airdata] = airdata_dynamic(alt, v)
    % computes static and dynamic air data from altitude and velocity vector, according to US standard atmosphere 
    % air data: density, dynamic pressure, mch number

    airdata = airdata_atmos(alt);
    airspeed = norm(v);

    % return values
    airdata.mach = airspeed / airdata.sonic_speed;
    airdata.dynamic_pressure = 0.5 * airdata.density * airspeed^2;
end