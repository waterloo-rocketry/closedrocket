function [meas_rotated] = imu_handler(meas, location)
    % Rotates IMU frame to body frame
    %#codegen
    
    persistent param
    if isempty(param)
        param = coder.load("model/model_params.mat");
    end

    meas_rotated = zeros(3,1);

    if location == "board"
        % rotate into body coordinates from sensor coordinate
        meas_rotated = param.S_board * meas;
    end
    if location == "mti"
        % rotate into body coordinates from sensor coordinate
        meas_rotated = param.S_mti * meas;
    end
    if location == "ad"
        % rotate into body coordinates from sensor coordinate
        meas_rotated = param.S_ad * meas;
    end
end