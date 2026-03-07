function [u, K] = controller_law_lqi(xR, r, L_delta, time)
    % tristan time :)
    % computes the optimal control signal for a flight condition 
    % u : control signal, desired canard angle (rad)
    % xR : roll state [roll angle (rad); roll rate (rad/s)]
    % r : reference signal, desired roll angle (rad)
    % L_delta : roll acceleration control derivative (rad/s^2 / rad)
    
    w_c = 1;

    persistent x_i r_f t
    if isempty(x_i)
        x_i = 0;
    end
    if isempty(r_f)
        r_f = r;
    end
    if isempty(t)
        t = -0.01;
    end

    r_f_dot = w_c * (r - r_f);
    r_f = r_f + (time - t) * (r_f_dot);
    x_i = x_i + (time - t) * (r_f - xR(1));
    t = time;
    z = [x_i; xR - [r_f; 0]]; % better performance than xR - [r_f; r_f_dot]
    


    %% augmented state, x = [int, phi, w] where int is the integral of roll error
    % A = [0  -1   0;
    %      0   0   1;
    %      0   0   0];
    % B = [0; 0; 1];
    % C = [0 1 0];
    % E = [1; 0; 0]; % integrator ff

    %% feedback gains
    % K = - 1 / L_delta * lqr(A, B, Q, R); % should test dlqr?
    % use lqr_calculator.m
    K = - 1 / L_delta * [-0.4472    3.7057    4.1727];

    %% control command
    u = K * z;
end
