function [x_new] = dynamics(dt, x, a)
    % Computes state update with dynamics model and time integration
    
    %% decomp
    % decompose state vector: [q(4); w(3); v(3); alt]
    q = x(1:4); w = x(5:7); v = x(8:10); alt = x(11);

    %% load parameters
    persistent param
    if isempty(param)
        param = coder.load("model_params.mat");
    end
    
    %% time updates
    % quaternion update
    q_new = quaternion_update(q, w, dt);

    % rate update
    w_exp_tilde = math_exp_tilde(w, dt);
    torque = dynamics_aero(v, alt, param);
    w_new = param.Jinv * (w_exp_tilde*param.J*w) + dt * param.Jinv * torque;
    % w_new = w + dt * param.Jinv * (torque - cross(w, param.J*w)); % old version
    %%% possibly more accurate: for Jx < Jy = Jz : u = (Jy-Jx)/Jy * wx, and 
    %%% wx_new = wx, [wy, wz]_new = Sx(u*dt)*[wy, wz] with Sx = [c(u), s(u); -s(u), c(u)]
    
    % velocity update 
    S = quaternion_rotmatrix(q); % inertial to body frame
    v_new = w_exp_tilde*v + dt * a + dt * S*param.g; 
    % v_new = v + dt * (a - cross(w,v) + S*param.g); % old version
        
    % altitude update
    v_earth = (S')*v;
    alt_new = alt + dt * v_earth(1);
    
    %% concoct state update vector
    x_new = [q_new; w_new; v_new; alt_new];
end