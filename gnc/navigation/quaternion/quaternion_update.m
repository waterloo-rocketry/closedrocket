function [q_new, dq] = quaternion_update(q_un, w, dt)
    % computes new quaternion from old quaternion and body rates

    %%% norm quaternions
    q = q_un / norm(q_un);
    
    %%% incremental quaternion
    dphi = norm(w) * dt / 2;
    if norm(w) == 0
        dn = [0; 0; 0];
    else    
        dn = w / norm(w);
    end
    dq = [cos(dphi); dn*sin(dphi)];
     
    %%% quaternion update
    q_new = quaternion_multiply(q, dq);
end