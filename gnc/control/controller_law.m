function [u, K, Kr] = controller_law(xR, r, L_delta)
    % computes the optimal control signal for a flight condition 
    % u : control signal, desired canard angle (rad)
    % xR : roll state [roll angle (rad); roll rate (rad/s)]
    % r : reference signal, desired roll angle (rad)
    % L_delta : roll acceleration control derivative (rad/s^2 / rad)
    
    %% tuning parameters
    %%% Q_phi : weight of angle error, Q_omega : weight of rate error
    %%% Mode 1: Tracking (low rates), Mode 2: Damping only (high rates)
    Q_phi_mode1   = 5;
    Q_omega_mode1 = 5;
    Q_phi_mode2   = 0;
    Q_omega_mode2 = 20;

    %%% thresholds
    w_low  = 0.5; % w < w_low: fully mode 1
    w_high = 1; % w > w_high: fully mode 2

    %% Control mode switching: linear crossfade
    w = abs(xR(2));
    blend = (w - w_low) / (w_high - w_low);
    blend = max(0, min(1, blend)); % clamp to [0,1]

    Q_phi = (1-blend)*Q_phi_mode1 + blend*Q_phi_mode2;
    Q_omega = (1-blend)*Q_omega_mode1 + blend*Q_omega_mode2;

    %% feedback gains
    K = - 1/L_delta * [ sqrt(Q_phi), sqrt( 2*sqrt(Q_phi) + Q_omega ) ];

    %% feedforward gain
    Kr = - K(1); % simplifies to this

    %% control command
    u = K * xR + Kr * r;
end
