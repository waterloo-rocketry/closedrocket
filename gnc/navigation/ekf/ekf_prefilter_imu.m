function [a, w] = ekf_prefilter_imu(board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, param)

    % base confidences (tune per sensor)

    % use accelerometer bias standard deviation
    C0_board_a = [1 1 1] / (1e-7)^2;
    C0_mti_a = [1 1 1] / (1e-7)^2;
    C0_ad_a = [1 1 1] / (1e-7)^2;

    % use gyroscope noise standard deviation
    C0_board_w = [1 1 1] / (1e-5)^2;
    C0_mti_w = [1 1 1] / (1e-5)^2;
    C0_ad_w = [1 1 1] / (1e-5)^2;

    % accel confidences (Hadamard product)
    C_board_a = C0_board_a .* board_accel.status;
    C_mti_a = C0_mti_a .* mti_accel.status;
    C_ad_a = C0_ad_a .* ad_accel.status;

    % gyro confidences (Hadamard product)
    C_board_w = C0_board_w .* board_gyro.status;
    C_mti_w = C0_mti_w .* mti_gyro.status;
    C_ad_w = C0_ad_w .* ad_gyro.status;

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

    % bias-corrected gyros
    w_board = board_gyro.gyro - board_gyro.bias;
    w_mti = mti_gyro.gyro - mti_gyro.bias;
    w_ad = ad_gyro.gyro - ad_gyro.bias;

    % weighted angular rate
    w = C_board_w .* w_board + C_mti_w .* w_mti + C_ad_w .* w_ad;

    % centrifugal correction
    a_board = board_accel.accel - cross(w, cross(w, param.d_board));
    a_mti = mti_accel.accel - cross(w, cross(w, param.d_mti));
    a_ad = ad_accel.accel - cross(w, cross(w, param.d_ad));

    % weighted acceleration
    a = C_board_a .* a_board + C_mti_a .* a_mti + C_ad_a .* a_ad;
end