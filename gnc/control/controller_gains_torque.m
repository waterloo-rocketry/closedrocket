function [K, Kr] = controller_gains_torque(L_delta, Q_phi, Q_omega)
    % computes the optimal control gains for a flight condition 
    % L_delta: roll acceleration control derivative, Q_phi: weight of angle error, Q_omega: weight of rate error
  
    %% feedback gains
    K = - 1/L_delta * [ sqrt(Q_phi), sqrt( 2*sqrt(Q_phi) + Q_omega ) ];

    %% feedforward gain
    % A = [0, 1; 0, 0];
    % B = [0; L_delta];
    % C = [1, 0]; % tracking channel
    % Kr = -1 / ( C * inv(A+B*K) * B );
    Kr = - K(1);
end