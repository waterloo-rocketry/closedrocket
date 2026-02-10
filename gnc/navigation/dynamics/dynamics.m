function [x_new] = dynamics(dt, x, a)
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
    torque = dynamics_aero(v, alt, param);
    
    %% time updates

    % quaternion update
    % q_new = quaternion_update(q, w, dt)
    q_new = quaternion_increment(q, w, dt);
    % q_err = q_new - quaternion_update(q, w, dt)

    % rate update
    w_new = w + dt * param.Jinv * (torque - cross(w, param.J*w));
    
    % velocity update 
    v_new = v + dt * (a - cross(w,v) + S*param.g);

    % altitude update
    v_earth = (S')*v;
    alt_new = alt + dt * v_earth(1);

    
    %% concoct state derivative vector
    x_new = [q_new; w_new; v_new; alt_new];
end
