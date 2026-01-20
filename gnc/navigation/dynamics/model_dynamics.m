function [x_new] = model_dynamics(dt, x, u)
    % Computes state derivative with predictive model. Use ODE solver to compute next state.
    
    %% decomp
    % decompose state vector: [q(4); w(3); v(3); alt]
    q = x(1:4); w = x(5:7); v = x(8:10); alt = x(11);

    % decompose input vector
    a = u.accel;

    %% load parameters
    persistent param
    if isempty(param)
        param = coder.load("model/model_params.mat");
    end
    
    %% compute rotation matrix 
    %%% attitude transformation, inertial to body frame
    S = quaternion_rotmatrix(q);

    %% forces and torques
    torque = aerodynamics(v, alt, param);
    
    %% time updates

    % quaternion update
    % q_new = quaternion_update(q, w, dt)
    q_new = quaternion_increment(q, w, dt);
    % q_err = q_new - quaternion_update(q, w, dt)

    % rate update
    w_new = w + dt * param.Jinv * (torque - cross(w, param.J*w));
    
    % velocity update 
    %%% acceleration specific force    
    v_new = v + dt * (a - cross(w,v) + S*param.g);

    % altitude update
    v_earth = (S')*v;
    alt_new = alt + dt * v_earth(1);

    
    %% concoct state derivative vector
    x_new = [q_new; w_new; v_new; alt_new];
end


%% aerodynamics
function [torque] = aerodynamics(v, alt, param)
    %%% air data  
    airdata = model_airdata(alt, v);
    p_dyn = airdata.dynamic_pressure;

    w = x(5:7); v = x(8:10);

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
