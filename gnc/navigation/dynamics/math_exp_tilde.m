function [skewed_exp] = math_exp_tilde(w, dt)
    % exp(-tilde(w)*dt) 
    dtw_tilde = math_tilde(-dt * w);
    skewed_exp = math_matrix_exp(dtw_tilde);
end