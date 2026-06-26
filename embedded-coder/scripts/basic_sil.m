% time = 12;
% dt_ctrl = 0.01;
% xR = [0.6; 0.4];
% pdyn = 80000;
% delta = 0.2;

ctrl_mem_in.coeffs = [2; 0];
ctrl_mem_in.w = 0.38;
ctrl_mem_in.P = diag([1e-9, 1e-9]);
ctrl_mem_in.delta_lp = 0.18;
ctrl_mem_in.w_dot_lp = 3.4;

[u, r, ctrl_mem_out, w_status_ctrl] = controller_codegen_entry(time, dt_ctrl, xR, pdyn, delta, ctrl_mem_in)
[u, r, ctrl_mem_out, w_status_ctrl] = GNC_codegen_SIL_sil('controller_codegen_entry', time, dt_ctrl, xR, pdyn, delta, ctrl_mem_in)

dt = 0.03;

flight_phase = true;

x = [
    0.826554201666454;
    0.297841698590636;
   -0.496402830984393;
   -0.099280566196879;
    3.5;
    0.7;
   -0.4;
    900;
    60;
    25;
    12000
];

P = [
    0.03, 0,    0,    0,    0,    0,    0,    0,    0,    0,    0;
    0,    0.03, 0,    0,    0,    0,    0,    0,    0,    0,    0;
    0,    0,    0.03, 0,    0,    0,    0,    0,    0,    0,    0;
    0,    0,    0,    0.03, 0,    0,    0,    0,    0,    0,    0;
    0,    0,    0,    0,    0.03, 0,    0,    0,    0,    0,    0;
    0,    0,    0,    0,    0,    0.03, 0,    0,    0,    0,    0;
    0,    0,    0,    0,    0,    0,    0.03, 0,    0,    0,    0;
    0,    0,    0,    0,    0,    0,    0,    0.03, 0,    0,    0;
    0,    0,    0,    0,    0,    0,    0,    0,    0.03, 0,    0;
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0.03, 0;
    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0.03
];

bias.board_gyro = [
    0.01;
    0.0;
   -0.01
];

bias.mti_gyro = [
    0.0;
    0.0;
    0.0
];

bias.ad_gyro = [
    0.0;
    0.0;
    0.0
];

bias.board_mag_earth = [
    0.2;
    0.0;
    0.44
];

bias.mti_mag_earth = [
    0;
    0;
    0
];

bias.board_baro = 99500.0;

bias.mti_baro = 0;

sens_filt.board_accel = [
    0.0;
    0.0;
    12.0
];

sens_filt.board_gyro = [
    0.0;
    0.0;
    0.0
];

sens_filt.mti_accel = [
    0;
    0;
    0
];

sens_filt.mti_gyro = [
    0;
    0;
    0
];

sens_filt.ad_accel = [
    0;
    0;
    0
];

sens_filt.ad_gyro = [
    0;
    0;
    0
];

sens_filt.board_baro = 99500.0;

sens_filt.board_mag = [
    0.2;
    0.0;
    0.44
];

sens_filt.mti_baro = 0;

sens_filt.mti_mag = [
    0;
    0;
    0
];

sens_input.board_accel.meas = [
    0.8;
   -0.4;
    48.0
];

sens_input.board_accel.status = true;

sens_input.board_gyro.meas = [
    0.35;
    0.12;
   -0.28
];

sens_input.board_gyro.status = true;

sens_input.mti_accel.meas = [
    0;
    0;
    0
];

sens_input.mti_accel.status = false;

sens_input.mti_gyro.meas = [
    0;
    0;
    0
];

sens_input.mti_gyro.status = false;

sens_input.ad_accel.meas = [
    0.6;
   -0.5;
    47.8
];

sens_input.ad_accel.status = true;

sens_input.ad_gyro.meas = [
    0.34;
    0.10;
   -0.27
];

sens_input.ad_gyro.status = true;

sens_input.board_baro.meas = 70000.0;

sens_input.board_baro.status = true;

sens_input.board_mag.meas = [
    0.16;
    0.08;
    0.39
];

sens_input.board_mag.status = true;

sens_input.mti_baro.meas = 0;

sens_input.mti_baro.status = false;

sens_input.mti_mag.meas = [
    0;
    0;
    0
];

sens_input.mti_mag.status = false;

[x_ret, P_ret, bias_ret, sens_filt_ret, cov_norm, roll_state, pdyn, w_status_nav] = navigation_codegen_entry(dt, flight_phase, x, P, bias, sens_filt, sens_input)

[x_ret, P_ret, bias_ret, sens_filt_ret, cov_norm, roll_state, pdyn, w_status_nav] = GNC_codegen_SIL_sil('navigation_codegen_entry', dt, flight_phase, x, P, bias, sens_filt, sens_input)
