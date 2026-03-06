function [R] = math_rotate(w, dt)
    % computes matrix exponential of rotation
    
    %%% incremental angle
    phi = norm(w) * dt;
    
    %%% normed skew-symmetric matrix
    if norm(w) == 0
        n = [0; 0; 0];
    else    
        n = w / norm(w);
    end
    n_tilde = math_tilde(n);

    %%% Rodrigues rotation formula
    R = eye(3) - sin(phi)*n_tilde + (1-cos(phi))*n_tilde^2;
end