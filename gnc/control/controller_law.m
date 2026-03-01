function [u, K, Kr] = controller_law(xR, r, L_delta)
    % computes the optimal control signal for a flight condition 
    % u : control signal, desired canard angle (rad)
    % xR : roll state [roll angle (rad); roll rate (rad/s)]
    % r : reference signal, desired roll angle (rad)
    % L_delta : roll acceleration control derivative (rad/s^2 / rad)

    %% control mode switching
    %%% mode = 1 : Tracking / low roll rate
    %%% mode = 2 : Damping / high roll rate
    persistent mode 
    if isempty(mode)
        mode = 1; % initialize
    end
    
    %%% switching hysteresis
    if mode == 1 && abs(xR(2)) >= 2 % upper threshold in rad/s
        mode = 2;
    elseif mode == 2 && ( ...
            abs(xR(2)) <= 0.5 && abs(xR(1)) <= 1 || ...
            abs(xR(2)) <= 0.2) % lower thresholds in rad/s and rad
        mode = 1;
    end
    
    %% tuning parameters
    %%% Q_phi : weight of angle error
    %%% Q_omega : weight of rate error
    switch mode
        case 1 % Control mode 1: Tracking / low roll rate
            Q_phi = 5; 
            Q_omega = 5;   
        otherwise % Control mode 2: Damping / high roll rate
            Q_phi = 0;
            Q_omega = 20;
    end 

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
