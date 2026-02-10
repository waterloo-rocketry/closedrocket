function [a, w] = ekf_prefilter_imu(dt, A1, W1, A2, W2, A3, W3)
    % Pre-filters redundant IMU data, weighted averages of acceleration and rates
    % Includes weighted averages of acceleration and rates, and bias handling
    % Inputs: time step dt, accelerometer Ai, gyroscope Wi 
    % Outputs: specific acceleration a, angular rates w

   %%% average specific force
    a = zeros(3,1);
    a1 = a;
    a2 = a;
    a_number = 0;
    if sensor_select(1) == 1 
        a1 = u(:,1) - cross(w, cross(w, param.d1)); % correction for centrifugal force
        a_number = a_number + 1;
    end
    if sensor_select(2) == 1
        a2 = u(:,2) - cross(w, cross(w, param.d2));
        a_number = a_number + 1;
    end
    if a_number ~= 0
        a = a1 + a2 / a_number; % average if multiple IMUs are alive
    end

    % decompose state vector: [q(4); w(3); v(3); alt]
    w = x(5:7); 
    % decompose bias matrix: [b_A(3,i); b_W(3, i); M_E(3, i); b_P(1, i)]
    b_W = b(4:6);
    y_expected = w + b_W; 

end