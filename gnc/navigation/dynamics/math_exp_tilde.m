function [skewed_exp] = math_exp_tilde(w, dt)
    % exp(-tilde(w)*dt) 
    skewed_exp = [1, exp(dt*w(3)), exp(-dt*w(2));
                  exp(-dt*w(3)), 1, exp(dt*w(1));
                  exp(dt*w(2)), exp(-dt*w(1)), 1];
end