function [qnew_q, qnew_w] = quaternion_update_jacobian(q, w, dt)
    % computes quaternion derivative from quaternion and body rates

    % norm quaternions
    q = quaternion_norm(q);

    % Quaternion product matrix
    Q = [q(1), -q(2), -q(3), -q(4);
         q(2), q(1), -q(4), q(3);
         q(3), q(4), q(1), -q(2);
         q(4), -q(3), q(2), q(1)];

    % Big Omega matrix
    W = [0, -w(1), -w(2), -w(3);
         w(1), 0, w(3), -w(2);
         w(2), -w(3), 0, w(1);
         w(3), w(2), -w(1), 0];

    % quaternion derivative Jacobians
    W_dt = 0.5 * dt * W;
    Q_dt =  0.5 * dt * Q;
    qnew_q = eye(4) + W_dt + 1/2 * W_dt^2 + 1/6 * W_dt^3 + 1/24 * W_dt^4;
    qnew_w = Q_dt;% + 1/2 * Q_dt^2 + 1/6 * Q_dt^3 + 1/24 * Q_dt^4;
    qnew_w = qnew_w(:, 2:4);
end