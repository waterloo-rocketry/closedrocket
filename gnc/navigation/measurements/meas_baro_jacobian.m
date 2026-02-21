function [J] = meas_baro_jacobian(x, bias)
    % Computes measurement prediction using current state and sensor biases

    %% decomp
    % decompose state vector: [q(4); w(3); v(3); alt;]
    alt = x(11);

    %% Initialize
    % Jacobian is of size: length(measurement) rows, length(x) columns
    % fill with zeros at first
    J = zeros(1, length(x));

    %% atmosphere model
    % [P, ~, ~, ~] = model_airdata(alt);
    airdata_altitude = airdata_atmos_jacobian(alt);
    P_alt = airdata_altitude.pressure;

    %% measurement prediction
    % y = [W; M; P];
    J(:, 11) = P_alt;
end