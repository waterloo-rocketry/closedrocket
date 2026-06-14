function [state, cov_norm, airdata, roll_state] = navigation_module(timestamp, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag)
    % Top-level navigation module. Calls code generation entry point.
    % Mocks firmware higher-level stuff

    persistent t flight_phase k first_run;

    %% Pad Filter
    %%% remember filtered values from last iteration
    persistent sens_filt x P bias

    %% initialize at beginning
    xhat = zeros(11,1);
    xhat(1) = 1;
    Phat = zeros(11);

    if isempty(x)
        x = xhat;
        P = Phat;
    end

    if isempty(sens_filt)
        sens_filt.board_accel = board_accel.meas;
        sens_filt.board_gyro = board_gyro.meas;
        sens_filt.mti_accel = mti_accel.meas;
        sens_filt.mti_gyro = mti_gyro.meas;
        sens_filt.ad_accel = ad_accel.meas;
        sens_filt.ad_gyro = ad_gyro.meas;
        sens_filt.board_baro = board_baro.meas;
        sens_filt.board_mag = board_mag.meas;
        sens_filt.mti_baro = mti_baro.meas;
        sens_filt.mti_mag = mti_mag.meas;
        bias.board_gyro = zeros(3,1);
        bias.mti_gyro = zeros(3,1);
        bias.ad_gyro = zeros(3,1);
        bias.board_mag_earth = zeros(3,1);
        bias.mti_mag_earth = zeros(3,1);
        bias.board_baro = 0;
        bias.mti_baro = 0;
    end

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
            flight_phase = true;
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

    sens_input.board_accel = board_accel;
    sens_input.board_gyro = board_gyro;
    sens_input.mti_accel = mti_accel;
    sens_input.mti_gyro = mti_gyro;
    sens_input.ad_accel = ad_accel;
    sens_input.ad_gyro = ad_gyro;
    sens_input.board_baro = board_baro;
    sens_input.board_mag = board_mag;
    sens_input.mti_baro = mti_baro;
    sens_input.mti_mag = mti_mag;

    [x, P, bias, sens_filt, airdata, roll_state] = navigation_codegen_entry(dt, flight_phase, x, P, bias, sens_filt, sens_input);

    %% Compute variance norm
    %%% for evaluating EKF quality
    cov_norm = norm(P); % scalar, 2-norm of the covariance matrix

    %% Pack state as struct
    %%% use union in C or smth
    state.q = x(1:4); % [1], attitude quaternion
    state.w = x(5:7); % [rad/s], angular rate
    state.v = x(8:10); % [m/s], velocity
    state.alt = x(11); % [m], altitude
    state.x = x; % also full state as vector is needed in simulink
end
