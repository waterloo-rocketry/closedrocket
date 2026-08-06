function [torque] = dynamics_aero(w, v, alt, param)
    % aerodynamics model

    %%% air data  
    airdata = airdata_atmos(alt);
    rho = airdata.density;

    %%% torques
    % signed quadratic cross flow, finite at zero forward speed and stays stable for negative forward vel
    torque_forcing = 0.5 * rho * param.c_aero * param.Cn_alpha ...
                   * [0; v(1)*v(3); -v(1)*v(2)];

    airspeed = norm(v);
    torque_damping = -0.5 * rho * airspeed ...
                   * param.c_aero_damping * param.Cn_omega ...
                   * [0; w(2); w(3)];

    torque = torque_forcing + torque_damping;
end
