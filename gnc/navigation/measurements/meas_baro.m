function [y] = meas_baro(x, bias)
    % Computes measurement prediction using current state and sensor biases

    %% decomp
    % decompose state vector: [q(4); w(3); v(3); alt]
    alt = x(11);

    % decompose bias matrix: [b_A(3,i); b_W(3, i); M_E(3, i); b_P(1, i)]
    b_P = bias(10);

    %% atmosphere model
    airdata = airdata_atmos(alt);
    P = airdata.pressure + b_P;

    %% measurement prediction
    y = [P];
end