function [u, K, Kr] = controller_law(xR, r, L_delta)
    % computes the optimal control signal for a flight condition 
    % u : control signal, desired canard angle (rad)
    % xR : roll state [roll angle (rad); roll rate (rad/s)]
    % r : reference signal, desired roll angle (rad)
    % L_delta : roll acceleration control derivative (rad/s^2 / rad)

    %% tuning parameters
    Q_phi = 10; % weight of angle error
    Q_omega = 10; % weight of rate error

    %% feedback gains
    K = - 1/L_delta * [ sqrt(Q_phi), sqrt( 2*sqrt(Q_phi) + Q_omega ) ];

    %% feedforward gain
    % A = [0, 1; 0, 0];
    % B = [0; L_delta];
    % C = [1, 0]; % tracking channel
    % Kr = -1 / ( C * inv(A+B*K) * B );
    Kr = - K(1); % simplifies to this

    %% control command
    u = K * xR + Kr * r;
end
