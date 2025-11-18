%% load
load("analysis\aurora\aurora_flight.mat");

%% plot est

%%% Colors (same as plot_state)
colors(1,:) = [0, 0, 0];      % Black
colors(2,:) = [0.9, 0.3, 0.1];% Red
colors(3,:) = [0.2, 0.8, 0.1];% Green
colors(4,:) = [0.1, 0.2, 0.8];% Blue

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')

f_q = figure(1);
plot(seconds(T.time_s), T.ATT_Q0, '.-', 'Color', colors(1,:), 'DisplayName', '$q_w$'); hold on;
plot(seconds(T.time_s), T.ATT_Q1, '.-', 'Color', colors(2,:), 'DisplayName', '$q_x$');
plot(seconds(T.time_s), T.ATT_Q2, '.-', 'Color', colors(3,:), 'DisplayName', '$q_y$');
plot(seconds(T.time_s), T.ATT_Q3, '.-', 'Color', colors(4,:), 'DisplayName', '$q_z$'); 
plot(seconds(TL.time_s), TL.ATT_Q0, 'Color', colors(1,:), 'HandleVisibility', 'off'); hold on;
plot(seconds(TL.time_s), TL.ATT_Q1, 'Color', colors(2,:), 'HandleVisibility', 'off');
plot(seconds(TL.time_s), TL.ATT_Q2, 'Color', colors(3,:), 'HandleVisibility', 'off');
plot(seconds(TL.time_s), TL.ATT_Q3, 'Color', colors(4,:), 'HandleVisibility', 'off');
% plot(seconds(TL.time_s), sqrt(TL.ATT_Q0.^2+TL.ATT_Q1.^2+TL.ATT_Q2.^2+TL.ATT_Q3.^2), 'k','HandleVisibility', 'off');
xlabel("Time [s]")
xlim([0, 140])
ylabel("Quaternion")
ylim([-1, 1])
% title("Attitude quaternion") 
legend('Location','southeast'); hold off;
grid on

f_e = figure(2);
plot(seconds(TL.time_s), TL.euler_roll, 'Color', colors(2,:), 'DisplayName', '$\phi$'); hold on;
plot(seconds(TL.time_s), TL.euler_pitch, 'Color', colors(3,:), 'DisplayName', '$\theta$');
plot(seconds(TL.time_s), TL.euler_yaw, 'Color', colors(4,:), 'DisplayName', '$\psi$');
xlabel("Time [s]")
xlim([0, 140])
ylabel("Angle [rad]")
% title("Relative Euler angles")
legend('Location','southwest'); hold off;
grid on

f_w = figure(3);
plot(seconds(T.time_s), T.RATE_WX, '.-', 'Color', colors(2,:), 'DisplayName', '$\omega_x$'); hold on;
plot(seconds(T.time_s), T.RATE_WY, '.-', 'Color', colors(3,:), 'DisplayName', '$\omega_y$')
plot(seconds(T.time_s), T.RATE_WZ, '.-', 'Color', colors(4,:), 'DisplayName', '$\omega_z$') 
plot(seconds(TL.time_s), TL.RATE_WX, 'Color', colors(2,:),'HandleVisibility', 'off');
plot(seconds(TL.time_s), TL.RATE_WY, 'Color', colors(3,:),'HandleVisibility', 'off');
plot(seconds(TL.time_s), TL.RATE_WZ, 'Color', colors(4,:),'HandleVisibility', 'off');
xlabel("Time [s]")
xlim([0, 140])
ylabel("Anglular velocity [rad/s]")
% title("Angular rates")
legend('Location','northwest'); hold off;
grid on

f_v = figure(4);
plot(seconds(T.time_s), T.VEL_VX, '.-', 'Color', colors(2,:),'DisplayName', '$v_x$'); hold on;
plot(seconds(T.time_s), T.VEL_VY, '.-', 'Color', colors(3,:),'DisplayName', '$v_y$');
plot(seconds(T.time_s), T.VEL_VZ, '.-', 'Color', colors(4,:),'DisplayName', '$v_z$');
plot(seconds(TL.time_s), TL.VEL_VX, 'Color', colors(2,:),'HandleVisibility', 'off');
plot(seconds(TL.time_s), TL.VEL_VY, 'Color', colors(3,:),'HandleVisibility', 'off');
plot(seconds(TL.time_s), TL.VEL_VZ, 'Color', colors(4,:),'HandleVisibility', 'off');
xlabel("Time [s]")
xlim([0, 140])
ylabel("Velocity [m/s]")
% title("Velocity")
legend('Location','best'); hold off;
grid on

f_a = figure(5);
plot(seconds(T.time_s), T.ALT/1000, '.-', 'Color', colors(4,:), 'DisplayName', '$r_x$'); hold on;
plot(seconds(TL.time_s), TL.ALT/1000, 'Color', colors(4,:),'HandleVisibility', 'off');
xlabel("Time [s]")
xlim([0, 140])
ylabel("Altitude [km]")
% title("Altitude")
legend('Location','best'); 
hold off;
grid on

f_c = figure(6);
plot(seconds(T.time_s), rad2deg(T.CANARD_ANGLE), '.-', 'Color', colors(2,:), 'DisplayName', '$\delta$'); hold on;
plot(seconds(TL.time_s),  rad2deg(TL.CANARD_ANGLE), 'Color', colors(2,:),'HandleVisibility', 'off');
plot(seconds(T.time_s), T.COEFF_CL, '.-', 'Color', colors(4,:), 'DisplayName', '$C_L$'), hold on;
plot(seconds(TL.time_s), TL.COEFF_CL,'Color', colors(4,:), 'HandleVisibility', 'off');
xlabel("Time [s]")
xlim([0, 140])
ylabel("Angle [deg], Coefficient")
% ylim([-1,5])
% title("Canard")
legend('Location','best'); hold off;
grid on

%%% plot control

f_cmd = figure(7);
plot(seconds(T_cmd.time_s), T_cmd.data, 'b-', 'DisplayName', 'cmd');
xlabel("Time [s]")
xlim([0, 140])
ylabel("Command [deg]")
ylim([-12,12])
% title("Canard")
%legend('Location','best'); hold off;
grid on

f_enc = figure(8);
plot(seconds(T_enc.time_s), T_enc.data, 'b.-', 'DisplayName', 'enc')
xlabel("Time [s]")
xlim([0, 140])
ylabel("Encoder [deg]")
ylim([-12,12])
% title("Canard")
%legend('Location','best'); hold off;
grid on


%%% plot IMU

f_imu_rate = figure(9);
plot(seconds(T_imu.time_s), T_imu.vel_x, 'r.-', 'DisplayName', '$\omega_x$'); hold on;
plot(seconds(T_imu.time_s), T_imu.vel_y, 'g.-', 'DisplayName', '$\omega_y$'); hold on;
plot(seconds(T_imu.time_s), T_imu.vel_z, 'b.-', 'DisplayName', '$\omega_z$'); hold on;
plot(seconds(TL_imu.time_s), TL_imu.vel_x, 'r','HandleVisibility', 'off');
plot(seconds(TL_imu.time_s), TL_imu.vel_y, 'g','HandleVisibility', 'off');
plot(seconds(TL_imu.time_s), TL_imu.vel_z, 'b','HandleVisibility', 'off');
xlabel("Time [s]")
xlim([0, 140])
ylabel("Anglular velocity [rad/s]")
legend('Location','best'); hold off;
grid on

f_imu_accel = figure(10);
plot(seconds(T_imu.time_s), T_imu.accel_x, 'r.-', 'DisplayName', '$a_x$'); hold on;
plot(seconds(T_imu.time_s), T_imu.accel_y, 'g.-', 'DisplayName', '$a_y$'); hold on;
plot(seconds(T_imu.time_s), T_imu.accel_z, 'b.-', 'DisplayName', '$a_z$'); hold on;
plot(seconds(TL_imu.time_s), TL_imu.accel_x, 'r','HandleVisibility', 'off');
plot(seconds(TL_imu.time_s), TL_imu.accel_y, 'g','HandleVisibility', 'off');
plot(seconds(TL_imu.time_s), TL_imu.accel_z, 'b','HandleVisibility', 'off');
xlabel("Time [s]")
xlim([0, 140])
ylabel("Acceleration [m/s$^2$]")
legend('Location','best'); hold off;
grid on


%%% plot rates
% f_rates = figure(11);
% plot(seconds(T_imu.time_s), T_imu.vel_x, 'r.-', 'DisplayName', 'x'); hold on;
% plot(seconds(T_imu.time_s), T_imu.vel_y, 'g.-', 'DisplayName', 'y')
% plot(seconds(T_imu.time_s), T_imu.vel_z, 'b.-', 'DisplayName', 'z')
% plot(seconds(T.time_s), T.RATE_WX, 'r.-', 'DisplayName', 'x');
% plot(seconds(T.time_s), T.RATE_WY, 'g.-', 'DisplayName', 'y');
% plot(seconds(T.time_s), T.RATE_WZ, 'b.-', 'DisplayName', 'z');
% xlabel("Time [s]")
% ylabel("Anglular rate [rad/s]")
% legend('Location','best'); hold off;
% grid on

set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')

fontsize(f_q, 12, "points")
fontsize(f_e, 12, "points")
fontsize(f_w, 12, "points")
fontsize(f_v, 12, "points")
fontsize(f_a, 12, "points")
fontsize(f_c, 12, "points")
fontsize(f_cmd, 12, "points")
fontsize(f_enc, 12, "points")
fontsize(f_imu_rate, 12, "points")
fontsize(f_imu_accel, 12, "points")


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
% exportgraphics(f_rates, 'analysis/aurora/aurora_rates.pdf')