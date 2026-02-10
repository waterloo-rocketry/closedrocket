function [torque_v] = dynamics_aero_jacobian(v, alt, param)
    % aerodynamics

    %%% air data 
    airdata = airdata_dynamic(alt, v);

    %torque_vx = Cl * delta * param.c_canard * [v(1), v(2), v(3); 0, 0, 0; 0, 0, 0];
    torque_vyz = 0.5 * param.c_aero * param.Cn_alpha * [0, 0, 0;
                                                        v(3), 0, v(1);
                                                        -v(2), -v(1), 0];
    %torque_v =  airdata.density * (torque_vx + torque_vyz);
    torque_v =  airdata.density * (torque_vyz);

end
