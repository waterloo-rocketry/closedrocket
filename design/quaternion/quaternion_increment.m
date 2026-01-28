function [q_new] = quaternion_increment(q_un, omega, dt)
    % computes new quaternion from old quaternion and body rates

    % norm quaternions
    q = q_un / norm(q_un);
    
    % incremental quaternion
    dphi = norm(omega) * dt / 2;
    if norm(omega) == 0
        dn = [0; 0; 0];
    else    
        dn = omega / norm(omega);
    end

    dq = [cos(dphi); dn*sin(dphi)];
     
    % quaternion derivative
    q_new = quaternion_multiply(q, dq);
end