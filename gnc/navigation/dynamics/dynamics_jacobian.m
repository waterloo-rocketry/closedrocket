function [J_x] = dynamics_jacobian(dt, x, a)
    % Computes state derivative with predictive model. Use ODE solver to compute next state.
    
    %% decomp
    % decompose state vector: [q(4); w(3); v(3)]
    q = x(1:4); w = x(5:7); v = x(8:10); alt = x(11);

    %% load parameters
    persistent param
    if isempty(param)
        param = coder.load("model_params.mat");
    end
    
    %% create empty Jacobian 
    J_x = zeros(length(x), length(x));
    % could also initialize as identity eye(length(x)), as right now all
    % sub-Jacobians have identity on the main diagonal
    
    %% quaternion attitude rows (q, 1:4)
    [q_q, q_w] = quaternion_update_jacobian(q, w, dt);

    J_x(1:4,1:4) = q_q; % column q (attitude)
    J_x(1:4, 5:7) = q_w; % column w (rates)

    %% angular rate rows (w, 5:7)
    % when implementing in C: the torque partial derivatives can probably be put in one function
    [torque_v] = dynamics_aero_jacobian(v, alt, param);

    w_w = eye(3) + dt * param.Jinv * tilde(param.J*w); % torque_w = 0 for now
    w_v = dt * param.Jinv * torque_v;

    J_x(5:7,5:7) = w_w; % column w
    J_x(5:7,8:10) = w_v; % column v


    %% velocity rows (v, 8:10)
    v_q = dt * quaternion_rotate_jacobian(q, param.g);
    v_w = dt * tilde(v);
    v_v = eye(3) - dt * tilde(w);

    J_x(8:10,1:4) = v_q; % column q
    J_x(8:10,5:7) = v_w; % column w
    J_x(8:10,8:10) = v_v; % column v


    %% altitude row (alt, 11)
    r_q = dt * quaternion_rotate_jacobian(quaternion_inv(q), v);
    alt_q = r_q(1,:); % only use altitude from position vector
    r_v = dt * quaternion_rotmatrix(quaternion_inv(q));
    alt_v = r_v(1,:);
    % r_r = eye(3);
    alt_alt = 1;

    J_x(11,1:4) = alt_q; % column q
    J_x(11,8:10) = alt_v; % column v
    J_x(11, 11) = alt_alt; % column alt

end