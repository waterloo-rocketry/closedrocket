function [a, w] = ekf_prefilter_imu(bias, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro)
    %%% computes average acceleration and angular rates from multiple IMUs.
    %%% includes correction of gyroscope bias and centrifugal acceleration.
    
    %% confidences
    %%% base confidences (tune per sensor)
    % use accelerometer bias standard deviation
    C0_board_a = [1; 1; 1] / (1e-7)^2;
    C0_mti_a = [1; 1; 1] / (1e-7)^2;
    C0_ad_a = [1; 1; 1] / (1e-7)^2;

    % use gyroscope noise standard deviation
    C0_board_w = [1; 1; 1] / (1e-5)^2;
    C0_mti_w = [1; 1; 1] / (1e-5)^2;
    C0_ad_w = [1; 0; 0] / (1e-5)^2;

    %%% confidence calculations
    % sensor status
    C_board_a = C0_board_a * board_accel.status;
    C_mti_a = C0_mti_a * mti_accel.status;
    C_ad_a = C0_ad_a * ad_accel.status;
    C_board_w = C0_board_w * board_gyro.status;
    C_mti_w = C0_mti_w * mti_gyro.status;
    C_ad_w = C0_ad_w * ad_gyro.status;

    C_total_a = C_board_a + C_mti_a + C_ad_a;
    C_total_w = C_board_w + C_mti_w + C_ad_w;

    if any(C_total_a == 0)
        error('No accel confidence on at least one axis.');
    end
    if any(C_total_w == 0)
        error('No gyro confidence on at least one axis.');
    end

    % normalize (Hadamard division)
    C_board_a = C_board_a ./ C_total_a;
    C_mti_a = C_mti_a ./ C_total_a;
    C_ad_a = C_ad_a ./ C_total_a;
    C_board_w = C_board_w ./ C_total_w;
    C_mti_w = C_mti_w ./ C_total_w;
    C_ad_w = C_ad_w ./ C_total_w;

    %% parameters
    persistent param
    if isempty(param)
        param = coder.load("gnc/model_params.mat");
    end 

    %% averaging
    % gyro bias correction
    w_board = board_gyro.meas - bias.board_gyro;
    w_mti = mti_gyro.meas - bias.mti_gyro;
    w_ad = ad_gyro.meas - bias.ad_gyro;

    % weighted angular rates
    w = C_board_w .* w_board + C_mti_w .* w_mti + C_ad_w .* w_ad;

    % centrifugal acceleration correction
    w_tilde = math_tilde(w);
    w_tilde_sq = w_tilde * w_tilde;
    a_board = board_accel.meas;% - w_tilde_sq * param.d_board;
    a_mti = mti_accel.meas;% - w_tilde_sq * param.d_mti;
    a_ad = ad_accel.meas;% - w_tilde_sq * param.d_ad;

    % weighted acceleration
    a = C_board_a .* a_board + C_mti_a .* a_mti + C_ad_a .* a_ad;
end