function [qnew_q, qnew_w] = quaternion_update_jacobian(q, w, dt)
    % computes quaternion derivative from quaternion and body rates

    q_norm = norm(q);
    q_unit = q / q_norm;

    angle = norm(w) * dt / 2;
    if norm(w) > 0
        dq = [cos(angle); w/norm(w)*sin(angle)];
    else
        dq = [1; 0; 0; 0];
    end

    % product matrix q*dq = R(dq)*q
    R_dq = [dq(1), -dq(2), -dq(3), -dq(4);
            dq(2),  dq(1),  dq(4), -dq(3);
            dq(3), -dq(4),  dq(1),  dq(2);
            dq(4),  dq(3), -dq(2),  dq(1)];
    norm_q = (eye(4) - q_unit*q_unit') / q_norm; %derivative of norm to remove radial mode
    qnew_q = R_dq * norm_q;

    % first order rate derivative
    Q = [q_unit(1), -q_unit(2), -q_unit(3), -q_unit(4);
         q_unit(2), q_unit(1), -q_unit(4), q_unit(3);
         q_unit(3), q_unit(4), q_unit(1), -q_unit(2);
         q_unit(4), -q_unit(3), q_unit(2), q_unit(1)];
    qnew_w = 0.5 * dt * Q(:,2:4);
end
