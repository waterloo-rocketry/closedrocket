function [x_new] = dynamics(dt, x, u)
    % Computes state derivative with predictive model. Use ODE solver to compute next state.
    
    %% decomp
    % decompose state vector: [q(4); w(3); v(3); alt]
    q = x(1:4); w = x(5:7); v = x(8:10); alt = x(11);

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
    %%% average specific force
    a = zeros(3,1);
    a1 = a;
    a2 = a;
    a_number = 0;
    if sensor_select(1) == 1 
        a1 = u(:,1) - cross(w, cross(w, param.d1)); % correction for centrifugal force
        a_number = a_number + 1;
    end
    if sensor_select(2) == 1
        a2 = u(:,2) - cross(w, cross(w, param.d2));
        a_number = a_number + 1;
    end
    if a_number ~= 0
        a = a1 + a2 / a_number; % average if multiple IMUs are alive
    end
    %%% acceleration  
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
