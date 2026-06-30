%% Controller single-point equivalence test

time = 12;
dt_ctrl = 0.01;
xR = [0.6; 0.4];          % [roll angle; roll rate]
pdyn_ctrl = 80000;        % Pa
delta = 0.2;              % rad

ctrl_mem_in.coeffs = [2; 0];
ctrl_mem_in.w = 0.38;
ctrl_mem_in.P = diag([1e-9, 1e-9]);
ctrl_mem_in.delta_lp = 0.18;
ctrl_mem_in.w_dot_lp = 3.4;

[u_m, r_m, ctrl_mem_m, w_status_ctrl_m] = ...
    controller_codegen_entry(time, dt_ctrl, xR, pdyn_ctrl, delta, ctrl_mem_in);

[u_s, r_s, ctrl_mem_s, w_status_ctrl_s] = ...
    GNC_codegen_SIL_sil('controller_codegen_entry', time, dt_ctrl, xR, pdyn_ctrl, delta, ctrl_mem_in);


%% Navigation single-point equivalence test

dt = 0.03;
flight_phase = true;

x = [
    0.815150960412541;
    0.293732638667334;
   -0.489554397778890;
   -0.097910879555778;
    3.5;
    0.7;
   -0.4;
    900;
    60;
    25;
    12000
];

P = eye(11) * 0.03;

bias.board_gyro = [0.01; 0.0; -0.01];
bias.mti_gyro = [0.0; 0.0; 0.0];
bias.ad_gyro = [0.0; 0.0; 0.0];

bias.board_mag_earth = [0.2; 0.0; 0.44];
bias.mti_mag_earth = [0; 0; 0];

bias.board_baro = 99500.0;
bias.mti_baro = 0;

sens_filt.board_accel = [0.0; 0.0; 12.0];
sens_filt.board_gyro = [0.0; 0.0; 0.0];

sens_filt.mti_accel = [0; 0; 0];
sens_filt.mti_gyro = [0; 0; 0];

sens_filt.ad_accel = [0; 0; 0];
sens_filt.ad_gyro = [0; 0; 0];

sens_filt.board_baro = 99500.0;
sens_filt.board_mag = [0.2; 0.0; 0.44];

sens_filt.mti_baro = 0;
sens_filt.mti_mag = [0; 0; 0];

sens_input.board_accel.meas = [0.8; -0.4; 48.0];
sens_input.board_accel.status = true;

sens_input.board_gyro.meas = [0.35; 0.12; -0.28];
sens_input.board_gyro.status = true;

sens_input.mti_accel.meas = [0; 0; 0];
sens_input.mti_accel.status = false;

sens_input.mti_gyro.meas = [0; 0; 0];
sens_input.mti_gyro.status = false;

sens_input.ad_accel.meas = [0.6; -0.5; 47.8];
sens_input.ad_accel.status = true;

sens_input.ad_gyro.meas = [0.34; 0.10; -0.27];
sens_input.ad_gyro.status = true;

sens_input.board_baro.meas = 70000.0;
sens_input.board_baro.status = true;

sens_input.board_mag.meas = [0.16; 0.08; 0.39];
sens_input.board_mag.status = true;

sens_input.mti_baro.meas = 0;
sens_input.mti_baro.status = false;

sens_input.mti_mag.meas = [0; 0; 0];
sens_input.mti_mag.status = false;

[x_m, P_m, bias_m, sens_filt_m, cov_norm_m, roll_state_m, pdyn_nav_m, w_status_nav_m] = ...
    navigation_codegen_entry(dt, flight_phase, x, P, bias, sens_filt, sens_input);

[x_s, P_s, bias_s, sens_filt_s, cov_norm_s, roll_state_s, pdyn_nav_s, w_status_nav_s] = ...
    GNC_codegen_SIL_sil('navigation_codegen_entry', dt, flight_phase, x, P, bias, sens_filt, sens_input);


%% Comparisons

absTol = 1e-10;
relTol = 1e-9;

compare_value('controller.u', u_m, u_s, absTol, relTol);
compare_value('controller.r', r_m, r_s, absTol, relTol);
compare_value('controller.ctrl_mem', ctrl_mem_m, ctrl_mem_s, absTol, relTol);
compare_value('controller.w_status_ctrl', w_status_ctrl_m, w_status_ctrl_s, absTol, relTol);

compare_value('navigation.x', x_m, x_s, absTol, relTol);
compare_value('navigation.P', P_m, P_s, absTol, relTol);
compare_value('navigation.bias', bias_m, bias_s, absTol, relTol);
compare_value('navigation.sens_filt', sens_filt_m, sens_filt_s, absTol, relTol);
compare_value('navigation.cov_norm', cov_norm_m, cov_norm_s, absTol, relTol);
compare_value('navigation.roll_state', roll_state_m, roll_state_s, absTol, relTol);
compare_value('navigation.pdyn', pdyn_nav_m, pdyn_nav_s, absTol, relTol);
compare_value('navigation.w_status_nav', w_status_nav_m, w_status_nav_s, absTol, relTol);

fprintf('\nAll MATLAB vs SIL comparisons passed.\n');


%% Local comparison helper

function compare_value(name, a, b, absTol, relTol)

    if isstruct(a)
        fields_a = fieldnames(a);
        fields_b = fieldnames(b);

        assert(isequal(sort(fields_a), sort(fields_b)), ...
            'Field mismatch in %s', name);

        for i = 1:numel(fields_a)
            field = fields_a{i};
            compare_value(name + "." + field, a.(field), b.(field), absTol, relTol);
        end

    elseif isnumeric(a)
        assert(isequal(size(a), size(b)), ...
            'Size mismatch in %s', name);

        diff = abs(a - b);
        tol = absTol + relTol * max(abs(a), abs(b));

        if any(diff(:) > tol(:))
            max_err = max(diff(:));
            error('Numeric mismatch in %s. Max error = %.3e', name, max_err);
        end

    elseif islogical(a)
        assert(isequal(a, b), ...
            'Logical mismatch in %s', name);

    else
        assert(isequal(a, b), ...
            'Value mismatch in %s', name);
    end
end

clear GNC_codegen_SIL_sil;