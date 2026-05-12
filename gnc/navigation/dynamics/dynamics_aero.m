function [torque] = dynamics_aero(v, alt, param)
    % aerodynamics model

    %%% air data  
    airdata = airdata_atmos(alt);
    airdata = airdata_dynamic(airdata, v);
    p_dyn = airdata.dynamic_pressure;

    %%% angle of attack/sideslip
    sin_alpha = sin(atan2(v(3),v(1)));
    sin_beta =  - sin(atan2(v(2),v(1)));

    %%% torques
    %torque_canards = Cl *  delta * param.c_canard * p_dyn *[1;0;0];
    torque_aero = p_dyn * ( param.c_aero * param.Cn_alpha * [0; sin_alpha; sin_beta] ); 
            %+ param.Cn_omega*[0; w(2); w(3)] ) * param.c_aero; % commented
            % out because timeline
    %torque = torque_canards + torque_aero;
    torque = torque_aero;
end
