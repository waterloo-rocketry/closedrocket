function [skewed_exp_vector] = math_exp_tilde_jacobian(w, dt, v)
    % jacobian of {exp(-tilde(w)*dt)*v} wrt w
    skewed_exp_vector = [0,                   -dt*v(3)*exp(-dt*w(2)), dt*v(2)*exp(w(3));
                        dt*v(3)*exp(dt*w(1)),  0,                    -dt*v(1)*exp(-dt*w(3));
                       -dt*v(2)*exp(-dt*w(1)), dt*v(1)*exp(dt*w(2)), 0];
end