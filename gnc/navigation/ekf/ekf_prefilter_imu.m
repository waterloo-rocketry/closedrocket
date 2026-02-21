function [a, w] = ekf_prefilter_imu(board_imu, mti_imu, ad_imu)
    % Pre-filters redundant IMU data, weighted averages of acceleration and rates
    % Includes weighted averages of acceleration and rates, and bias handling
    % Inputs: time step dt, accelerometer Ai, gyroscope Wi 
    % Outputs: specific acceleration a, angular rates w
    % Hadamard product / division is simply performed element by element
    
    % sensor confidences
    C_board_a = [1 1 1] * 10;
    C_board_a = C_board_a .* board_imu.accel_status; % Hadamard product
    C_board_w = [1 1 1] * 10;
    C_board_w = C_board_w .* board_imu.gyro_status;

    C_mti_a = [1 1 1] * 5;
    C_mti_a = C_mti_a .* mti_imu.accel_status;
    C_mti_w = [1 1 1] * 5;
    C_mti_w = C_mti_w .* mti_imu.gyro_status;

    C_ad_a = [1 1 1] * 1;
    C_ad_a = C_ad_a .* ad_imu.accel_status;
    C_ad_w = [1 0 0] * 1;
    C_ad_w = C_ad_w .* ad_imu.gyro_status;

    C_total_a = C_board_a + C_mti_a + C_ad_a;
    C_total_w = C_board_w + C_mti_w + C_ad_w;

    % if zero confidence in all sensors for any parameter
    % should probably introduce actual handling (or start crying lol)
    if any(C_total_a == 0)
        error('No accel confidence on at least one axis.');
    end
    if any(C_total_w == 0)
        error('No gyro confidence on at least one axis.');
    end

    % normalize confidence with Hadamard division
    C_board_a = C_board_a ./ C_total_a;
    C_mti_a = C_mti_a ./ C_total_a;
    C_ad_a = C_ad_a ./ C_total_a;
    
    C_board_w = C_board_w ./ C_total_w;
    C_mti_w = C_mti_w ./ C_total_w;
    C_ad_w = C_ad_w ./ C_total_w;

    
    %%% angular rates bias correction
    w_board = board_imu.gyro - board_imu.gyro_bias;
    w_mti = mti_imu.gyro - mti_imu.gyro_bias;
    w_ad = ad_imu.gyro - ad_imu.gyro_bias;

    %%% angular rates average
    w = C_board_w .* w_board + C_mti_w .* w_mti + C_ad_w .* w_ad;

    %%% centrifugal correction
    a_board = board_imu.accel - cross(w, cross(w, param.d_board));
    a_mti = mti_imu.accel - cross(w, cross(w, param.d_mti));
    a_ad = ad_imu.accel - cross(w, cross(w, param.d_ad));

    %%% acceleration average
    a = C_board_a .* a_board + C_mti_a .* a_mti + C_ad_a .* a_ad;

end