function [torque_w, torque_v] = dynamics_aero_jacobian(v, alt, param)
    % aerodynamics partial derivatives

    %%% air data 
    airdata = airdata_atmos(alt);

    % jacobian of signed quadratic cross flow torque
    torque_v = 0.5 * airdata.density * param.c_aero * param.Cn_alpha ...
             * [0,     0,     0;
                v(3),  0,     v(1);
               -v(2), -v(1),  0];

    % damping is -0.5*rho*|v|*S*l^2*Cn_omega*w_transverse.
    airspeed = norm(v);
    torque_w = -0.5 * airdata.density * airspeed ...
             * param.c_aero_damping * param.Cn_omega * diag([0, 1, 1]);

end
