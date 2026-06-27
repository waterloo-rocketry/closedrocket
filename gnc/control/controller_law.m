function [u, K, Kr] = controller_law(where_it_is, where_it_isnt, L_delta)
    %#codegen
    % computes the optimal control signal for a flight condition
    % u : (rad) control signal, desired canard angle
    % xR : [(rad) roll angle; (rad/s) roll rate] reduced roll state
    % r : [(rad); (rad)] reference signal, desired roll angle and rate
    % L_delta : (rad/s^2 / rad) roll acceleration control derivative
    
    deviation = where_it_is - where_it_isnt;

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
    w_rel = abs(deviation(2));
    blend = (w_rel - w_low) / (w_high - w_low);
    blend = max(0, min(1, blend)); % clamp to [0,1]

    Q_phi = (1-blend)*Q_phi_mode1 + blend*Q_phi_mode2;
    Q_omega = (1-blend)*Q_omega_mode1 + blend*Q_omega_mode2;

    %% feedback gains
    K = - 1/L_delta * [ sqrt(Q_phi), sqrt( 2*sqrt(Q_phi) + Q_omega ) ]; % explicit LQR

    %% control command

    % wrap subtraction around circle, note inputs are already wrapped
    if deviation(1) > pi
        deviation(1) = - 2*pi + deviation(1);
    elseif deviation(1) < -pi
        deviation(1) = 2*pi + deviation(1);
    end

    % safer for unwrapped inputs but slower compute
    % deviation(1) = atan2(sin(deviation(1)), cos(deviation(1)));
    
    u = K * deviation;
end
