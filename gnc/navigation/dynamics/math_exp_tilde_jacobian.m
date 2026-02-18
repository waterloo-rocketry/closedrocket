function [skewed_exp_w] = math_exp_tilde_jacobian(w, dt, v)
    % jacobian of {exp(-tilde(w)*dt)*v} wrt w   
    skewed_exp_w = dt*math_tilde(v) +  1/2*dt^2 * (math_tilde(v)*math_tilde(w) - 2*math_tilde(w)*math_tilde(v));
end