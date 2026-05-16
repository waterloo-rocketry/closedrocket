function [y] = meas_baro(x, b_P)
    % Computes measurement prediction using current state and sensor biases

    %% decomp
    % decompose state vector: [q(4); w(3); v(3); alt]
    alt = x(11);

    %% atmosphere model
    airdata = airdata_atmos(alt);
    P = airdata.pressure + b_P; % [Pa], measured atmospheric pressure

    %% measurement prediction
    y = [P];
end