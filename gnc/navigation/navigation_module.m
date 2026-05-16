function [state, cov_norm, airdata, roll_state] = navigation_module(timestamp, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag)
    % Top-level navigation module. Calls code generation entry point.  
    % Mocks firmware higher-level stuff
    
    persistent t flight_phase k first_run;
    
    %% config settings
    idle_time = 15; % [s], wait time to handover from pad to flight, 5s before liftoff
    sampling_imu = 0.002; % [s], sampling period of imu
    sampling_other = 0.02; % [s], sampling period of baro, mag

    if isempty(first_run)
        first_run = 1;
        flight_phase = false;
        k = 1;
        if timestamp >= 0.005
                t = timestamp-0.005;
        else 
                t = 0;
        end
    end
    
    %% timecode
    dt = timestamp - t; % [s], time step size for dynamics integration
    t = timestamp;
    
    %% flight phase
    if t >= idle_time % mock for flight phase
            flight_phase = true; % 1 is pad, 0 is flight
    end

    %% Mock slower sampling rate of baro, mag
    if k ==  sampling_other/sampling_imu
        k = 1; % reset counter for the next cycle
    else
        k = k + 1; % increment counter for the next cycle
        board_baro.status = false; % sensors do not have new data
        board_mag.status = false;
        mti_baro.status = false;
        mti_mag.status = false;
    end

    [state, cov_norm, airdata, roll_state] = navigation_codegen_entry(dt, flight_phase, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag);
end