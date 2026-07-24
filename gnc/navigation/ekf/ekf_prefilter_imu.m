function [a, w, status_fast] = ekf_prefilter_imu(bias, sens_in)
    %%% computes average acceleration and angular rates from multiple IMUs.
    %%% includes correction of gyroscope bias and centrifugal acceleration.

    status_fast = false;
    a = zeros(3, 1);
    w = zeros(3, 1);
    
    %% parameters
    persistent param
    if isempty(param)
        param = coder.load("gnc/model_params.mat");
    end


    %% confidences
    %%% base confidences (tune per sensor)
    % use accelerometer bias standard deviation
    C0_board_a = [1; 1; 1] / (1e-7)^2;
    C0_mti_a = [1; 1; 1] / (1e-7)^2;
    C0_ad_a = [1; 1; 1] / (1e-7)^2;

    % use gyroscope noise standard deviation
    C0_board_w = [1; 1; 1] / (1e-5)^2;
    C0_mti_w = [1; 1; 1] / (1e-5)^2;
    C0_ad_w = [1; 0; 0] / (1e-3)^2;

    %% confidence calculations
    % sensor status
    C_board_a = C0_board_a * sens_in.board_accel.status;
    C_mti_a = C0_mti_a * sens_in.mti_accel.status;
    C_ad_a = C0_ad_a * sens_in.ad_accel.status;
    C_board_w = C0_board_w * sens_in.board_gyro.status;
    C_mti_w = C0_mti_w * sens_in.mti_gyro.status;
    C_ad_w = C0_ad_w * sens_in.ad_gyro.status;

    C_total_a = C_board_a + C_mti_a + C_ad_a;
    C_total_w = C_board_w + C_mti_w + C_ad_w;

    if any(C_total_a == 0) || any(C_total_w == 0)
        return;
    end

    status_fast = true;

    % normalize (Hadamard division)
    C_board_a = C_board_a ./ C_total_a;
    C_mti_a = C_mti_a ./ C_total_a;
    C_ad_a = C_ad_a ./ C_total_a;
    C_board_w = C_board_w ./ C_total_w;
    C_mti_w = C_mti_w ./ C_total_w;
    C_ad_w = C_ad_w ./ C_total_w;


    %% averaging
    % gyro bias correction
    w_board = sens_in.board_gyro.meas - bias.board_gyro;
    w_mti = sens_in.mti_gyro.meas - bias.mti_gyro;
    w_ad = sens_in.ad_gyro.meas - bias.ad_gyro;

    % weighted angular rates
    w = C_board_w .* w_board + C_mti_w .* w_mti + C_ad_w .* w_ad; % [rad/s]

    % centrifugal acceleration correction
    w_tilde = math_tilde(w);
    w_tilde_sq = w_tilde * w_tilde;
    a_board = sens_in.board_accel.meas - w_tilde_sq * param.d_board;
    a_mti = sens_in.mti_accel.meas - w_tilde_sq * param.d_mti;
    a_ad = sens_in.ad_accel.meas - w_tilde_sq * param.d_ad;

    % weighted acceleration
    a = C_board_a .* a_board + C_mti_a .* a_mti + C_ad_a .* a_ad; % [m/s^2]

end
