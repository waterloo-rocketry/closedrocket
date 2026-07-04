function [meas_rotated] = sensor_handler(meas, location, type)
    % Rotates IMU frame to body frame
    %#codegen
    
    persistent param
    if isempty(param)
        param = coder.load("gnc/model_params.mat");
    end

    meas_rotated = zeros(3,1);

    if location == "board"
        % rotate into body coordinates from sensor coordinate
        if type == "imu"
            meas_rotated = param.S_board_acc * meas;
        end
        if type == "mag"
            meas_rotated = param.S_board_mag * meas;
        end
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