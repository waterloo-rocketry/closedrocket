function [a, w] = ekf_prefilter_imu(dt, A1, W1, A2, W2, A3, W3)
    % Pre-filters redundant IMU data, weighted averages of acceleration and rates
    % Includes weighted averages of acceleration and rates, and bias handling
    % Inputs: time step dt, accelerometer Ai, gyroscope Wi 
    % Outputs: specific acceleration a, angular rates w

    a = zeros(3,1);
    w = zeros(3,1);

    %%% gyroscope bias correction
    % w1 = W1 - b_w1;

    %%% average angular rates
    % r1 = 1e-7; % use actual noise variances from datasheets
    % r2 = 1e-5; 
    % if sensor1_isdead == 1
    %     r1 = 1e10; % very very high
    % end
    % R1 = 1 - r1 / (r1+r2);
    %
    % w = R1* w1 + R2 * w2;

    %%% centrifugal correction
    % a1 = A1 - cross(w, cross(w, param.d1));
    
    %%% average acceleration
    % r1 = 1e-7; % use actual noise variances from datasheets
    % r2 = 1e-5; 
    % if sensor1_isdead == 1
    %     r1 = 1e10; % very very high
    % end
    % R1 = 1 - r1 / (r1+r2);
    %
    % a = R1*a1 + R2 *a2;

end