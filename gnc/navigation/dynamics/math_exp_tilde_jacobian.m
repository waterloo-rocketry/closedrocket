function [skewed_exp_v] = math_exp_tilde_jacobian(w, dt, v)
    % jacobian of {exp(-tilde(w)*dt)*v} wrt w
    dtw_tilde = math_tilde(-dt * w);
    skewed_exp = math_matrix_exp(dtw_tilde);
    jacobian_left = eye(3) - 1/2 * dtw_tilde + 1/6 * dtw_tilde^2 + 1/24*dtw_tilde^3;
    skewed_exp_v = - dt * skewed_exp * math_tilde(v) * jacobian_left;
end