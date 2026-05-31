function [state, cov_norm, airdata, roll_state] = navigation_module(timestamp, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag)
    % Top-level navigation module. Calls code generation entry point.  
    % Mocks firmware higher-level stuff
    
    persistent t flight_phase k first_run;

    %% Pad Filter
    %%% remember filtered values from last iteration
    persistent sf x P bias

    %% initialize at beginning
    xhat = zeros(11,1); 
    xhat(1) = 1; 
    Phat = zeros(11);

    if isempty(x)
        x = xhat; 
        P = Phat;
    end

    if isempty(sf)
        sf.board_accel_f = board_accel.meas;
        sf.board_gyro_f = board_gyro.meas;
        sf.mti_accel_f = mti_accel.meas;
        sf.mti_gyro_f = mti_gyro.meas; 
        sf.ad_accel_f = ad_accel.meas; 
        sf.ad_gyro_f = ad_gyro.meas;  
        sf.board_baro_f = board_baro.meas; 
        sf.board_mag_f = board_mag.meas; 
        sf.mti_baro_f = mti_baro.meas; 
        sf.mti_mag_f = mti_mag.meas; 
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

    [x_ret, P_ret, b_ret, sf_ret] = navigation_codegen_entry(dt, flight_phase, x, P, bias, sf, board_accel, board_gyro, mti_accel, mti_gyro, ad_accel, ad_gyro, board_baro, board_mag, mti_baro, mti_mag);

    x = x_ret;
    P = P_ret;
    bias = b_ret;

    %% Pack state as struct
    %%% use union in C or smth
    state.q = x(1:4); % [1], attitude quaternion
    state.w = x(5:7); % [rad/s], angular rate 
    state.v = x(8:10); % [m/s], velocity
    state.alt = x(11); % [m], altitude
    state.x = x; % also full state as vector is needed in simulink

    %% Compute variance norm 
    %%% for evaluating EKF quality
    cov_norm = norm(P); % scalar, 2-norm of the covariance matrix

    %% Compute air data
    airdata = airdata_atmos(x(11));
    airdata = airdata_dynamic(airdata, x(8:10));

    %% controller input vector
    phi = quaternion_to_roll(x(1:4));
    roll_state = [phi; x(5)];

    sf.board_accel_f = sf_ret.board_accel_f;
    sf.board_gyro_f = sf_ret.board_gyro_f;
    sf.mti_accel_f = sf_ret.mti_accel_f;
    sf.mti_gyro_f = sf_ret.mti_gyro_f;
    sf.ad_accel_f = sf_ret.ad_accel_f;
    sf.ad_gyro_f = sf_ret.ad_gyro_f;
    sf.board_baro_f = sf_ret.board_baro_f;
    sf.board_mag_f = sf_ret.board_mag_f;
    sf.mti_baro_f = sf_ret.mti_baro_f;
    sf.mti_mag_f = sf_ret.mti_mag_f;
    
end