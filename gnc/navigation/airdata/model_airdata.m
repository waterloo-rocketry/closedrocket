function [air_data] = model_airdata(alt, v)
    % computes static and dynamic air data from altitude and velocity vector, according to US standard atmosphere 
    % air data: density, dynamic pressure, mch number

    atmos_data = model_atmos(alt);
    airspeed = norm(v);

    % return values
    air_data.density = atmos_data.density;
    air_data.mach = airspeed / atmos_data.mach;
    air_data.dynamic_pressure = 0.5 * atmos_data.density * airspeed^2;
end