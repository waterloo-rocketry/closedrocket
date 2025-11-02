clear
% Read data (if not already in workspace)
% T_long = readtable('analysis/aurora/aurora_flight_data.xlsx'); 
T_est = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'state_est_data'); 
T_cmd = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'proc_cmd'); 
T_enc = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'mcb_encoder'); 
T_imu = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'altimu_meas'); 

%%

% relevant times only
% keep only times between time_start and time_end
time_start = 0;
time_end = 150;
time_proc_offset = 10863;
time_mcb_offset = 30644;% -2e3;

% Convert long to wide format
T = unstack(T_est, 'data', 'state_id');
% format names
oldNames = T.Properties.VariableNames;
newNames = strrep(oldNames, 'STATE_ID_', '');
T.Properties.VariableNames = newNames;

% replace NaN with previous 
TF = fillmissing(T, 'linear');

% Convert from milliseconds to seconds
function T = retimer(T, time_offset, time_start, time_end)
    T.time_s = (T.time_ms - time_offset)/ 1000;
    T.time_ms = [];
    % T = T(T.time_s >= time_start & T.time_s <= time_end, :);
    % T.time_s = T.time_s - T.time_s(1);
end
T = retimer(T, time_proc_offset, time_start, time_end);
TF = retimer(TF, time_proc_offset, time_start, time_end);
T_cmd = retimer(T_cmd, time_proc_offset, time_start, time_end);
T_enc = retimer(T_enc, time_mcb_offset, time_start, time_end);
T_imu = retimer(T_imu, time_proc_offset, time_start, time_end);


%% process data
%%% euler angles
for i=1:height(T)
    q = [T.ATT_Q0(i), T.ATT_Q1(i), T.ATT_Q2(i), T.ATT_Q3(i)]';
    euler = quaternion_to_euler(q);
    T.euler_roll(i) = euler(1);
    T.euler_pitch(i) = euler(2);
    T.euler_yaw(i) = euler(3);
end

%%% command, encoder
T_cmd.data = (T_cmd.data - 32768) / 1000;
T_enc.data = (T_enc.data - 32768) / 1000;


%%% IMU data 
S2 = [0, 0, -1;
     -1, 0, 0;
      0, 1, 0];
% to signed variables
signed = @(u) int16(bitset(u,16,0)) + (-2^15)*int16(bitget(u,16));
for i = 3:8
    for k = 1:22
        var = T_imu.(i)(k);
        if ~isnan(var)
            T_imu.(i)(k) = signed(var);
        end
    end
end
factor_gyro = deg2rad(2000 / 32768);
factor_accel = 16*9.81 / 32768;

T_imu.vel_x = T_imu.vel_x * factor_gyro;
T_imu.vel_y = T_imu.vel_y * factor_gyro;
T_imu.vel_z = T_imu.vel_z * factor_gyro;
T_imu.accel_x = T_imu.accel_x * factor_accel;
T_imu.accel_y = T_imu.accel_y * factor_accel;
T_imu.accel_z = T_imu.accel_z * factor_accel;

for k = 1:22
    rate = [-T_imu.vel_z(k); -T_imu.vel_x(k); T_imu.vel_y(k)];
    accel = [-T_imu.accel_z(k); -T_imu.accel_x(k); T_imu.accel_y(k)];
    T_imu.vel_x(k) = rate(1); 
    T_imu.vel_y(k) = rate(2); 
    T_imu.vel_z(k) = rate(3);
    T_imu.accel_x(k) = accel(1); 
    T_imu.accel_y(k) = accel(2); 
    T_imu.accel_z(k) = accel(3);
end

%%% dispersion and recovery
q_apo = [T.ATT_Q0(12); T.ATT_Q1(12); T.ATT_Q2(12); T.ATT_Q3(12)];
vel_apo_b = [T.VEL_VX(12); T.VEL_VY(12); T.VEL_VZ(12)];
vel_apo_g = quaternion_rotmatrix(q_apo) * vel_apo_b;
vel_apo_vert = vel_apo_g(1);
vel_apo_hor = norm(vel_apo_g(2:3));

%% plot est

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')

f_q = figure(1);
plot(T.time_s, T.ATT_Q0, 'o:', 'DisplayName', '$q_w$'); hold on;
plot(T.time_s, T.ATT_Q1, 'o:', 'DisplayName', '$q_x$');
plot(T.time_s, T.ATT_Q2, 'o:', 'DisplayName', '$q_y$');
plot(T.time_s, T.ATT_Q3, 'o:', 'DisplayName', '$q_z$');
% plot(TF.time_s, TF.ATT_Q0, 'o:', 'DisplayName', 'w'); hold on;
% plot(TF.time_s, TF.ATT_Q1, 'o:', 'DisplayName', 'x');
% plot(TF.time_s, TF.ATT_Q2, 'o:', 'DisplayName', 'y');
% plot(TF.time_s, TF.ATT_Q3, 'o:', 'DisplayName', 'z');
xlabel("Time [s]")
ylabel("Quaternion [ ]")
ylim([-1, 1])
% title("Attitude quaternion") 
legend('Location','southeast'); hold off;
grid on

f_e = figure(2);
plot(T.time_s, T.euler_roll, 'o:', 'DisplayName', '$\phi$'); hold on;
plot(T.time_s, T.euler_pitch, 'o:', 'DisplayName', '$\theta$');
plot(T.time_s, T.euler_yaw, 'o:', 'DisplayName', '$\psi$');
xlabel("Time [s]")
ylabel("Angle [rad]")
% title("Relative Euler angles")
legend('Location','southwest'); hold off;
grid on

f_w = figure(3);
plot(T.time_s, T.RATE_WX, 'o:', 'DisplayName', '$\omega_x$'); hold on;
plot(T.time_s, T.RATE_WY, 'o:', 'DisplayName', '$\omega_y$')
plot(T.time_s, T.RATE_WZ, 'o:', 'DisplayName', '$\omega_z$')
xlabel("Time [s]")
ylabel("Anglular rate [rad/s]")
% title("Angular rates")
legend('Location','northwest'); hold off;
grid on

f_v = figure(4);
plot(T.time_s, T.VEL_VX, 'o:', 'DisplayName', '$v_x$'); hold on;
plot(T.time_s, T.VEL_VY, 'o:', 'DisplayName', '$v_y$');
plot(T.time_s, T.VEL_VZ, 'o:', 'DisplayName', '$v_z$');
xlabel("Time [s]")
ylabel("Velocity [m/s]")
% title("Velocity")
legend('Location','best'); hold off;
grid on

f_a = figure(5);
plot(T.time_s, T.ALT, 'o:', 'DisplayName', '$r_x$')
xlabel("Time [s]")
ylabel("Altitude [m]")
% title("Altitude")
legend('Location','best'); 
hold off;
grid on

f_c = figure(6);
plot(T.time_s, rad2deg(T.CANARD_ANGLE), 'o:', 'DisplayName', '$\delta$'); hold on;
plot(T.time_s, T.COEFF_CL, 'o:', 'DisplayName', '$C_L$')
xlabel("Time [s]")
ylabel("Angle [deg], Coefficient [ ]")
% ylim([-1,5])
% title("Canard")
legend('Location','best'); hold off;
grid on

%% plot control

f_cmd = figure(7);
plot(T_cmd.time_s, T_cmd.data, '.:', 'DisplayName', 'cmd');
xlabel("Time [s]")
ylabel("Command [deg]")
ylim([-12,12])
% title("Canard")
%legend('Location','best'); hold off;
grid on

f_enc = figure(8);
plot(T_enc.time_s, T_enc.data, 'o-', 'DisplayName', 'enc')
xlabel("Time [s]")
ylabel("Encoder [deg]")
ylim([-12,12])
% title("Canard")
%legend('Location','best'); hold off;
grid on


%% plot IMU

f_imu_rate = figure(9);
plot(T_imu.time_s, T_imu.vel_x, 'o:', 'DisplayName', '$\omega_x$'); hold on;
plot(T_imu.time_s, T_imu.vel_y, 'o:', 'DisplayName', '$\omega_y$')
plot(T_imu.time_s, T_imu.vel_z, 'o:', 'DisplayName', '$\omega_z$')
xlabel("Time [s]")
ylabel("Anglular rate [rad/s]")
legend('Location','best'); hold off;
grid on

f_imu_accel = figure(10);
plot(T_imu.time_s, T_imu.accel_x, 'o:', 'DisplayName', '$a_x$'); hold on;
plot(T_imu.time_s, T_imu.accel_y, 'o:', 'DisplayName', '$a_y$')
plot(T_imu.time_s, T_imu.accel_z, 'o:', 'DisplayName', '$a_z$')
xlabel("Time [s]")
ylabel("Acceleration [m/s$^2$]")
legend('Location','best'); hold off;
grid on


%% plot rates
f_rates = figure(11);
plot(T_imu.time_s, T_imu.vel_x, 'ro:', 'DisplayName', 'x'); hold on;
plot(T_imu.time_s, T_imu.vel_y, 'go:', 'DisplayName', 'y')
plot(T_imu.time_s, T_imu.vel_z, 'bo:', 'DisplayName', 'z')
plot(T.time_s, T.RATE_WX, 'ro:', 'DisplayName', 'x');
plot(T.time_s, T.RATE_WY, 'go:', 'DisplayName', 'y');
plot(T.time_s, T.RATE_WZ, 'bo:', 'DisplayName', 'z');
xlabel("Time [s]")
ylabel("Anglular rate [rad/s]")
legend('Location','best'); hold off;
grid on

set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')


%% print Recovery
fprintf("Recovery: at %f s and %f m altitude \n the velocity was %f m/s vertical, %f m/s lateral.\n", T.time_s(12), T.ALT(12), vel_apo_vert, vel_apo_hor);

%% export
exportgraphics(f_q, 'analysis/aurora/aurora_q.pdf')
exportgraphics(f_e, 'analysis/aurora/aurora_euler.pdf')
exportgraphics(f_w, 'analysis/aurora/aurora_w.pdf')
exportgraphics(f_v, 'analysis/aurora/aurora_v.pdf')
exportgraphics(f_a, 'analysis/aurora/aurora_alt.pdf')
exportgraphics(f_c, 'analysis/aurora/aurora_canard.pdf')
exportgraphics(f_cmd, 'analysis/aurora/aurora_cmd.pdf')
exportgraphics(f_enc, 'analysis/aurora/aurora_enc.pdf')
exportgraphics(f_imu_rate, 'analysis/aurora/aurora_imu_rate.pdf')
exportgraphics(f_imu_accel, 'analysis/aurora/aurora_imu_accel.pdf')
exportgraphics(f_rates, 'analysis/aurora/aurora_rates.pdf')