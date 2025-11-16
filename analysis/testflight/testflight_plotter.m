%% load
load("analysis\testflight\testflight_flight.mat")

%% plot

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')

f_q = figure(1);
plot(T.timestamp_s, T.ATT_Q0, 'k', 'DisplayName', '$q_w$'); hold on;
plot(T.timestamp_s, T.ATT_Q1, 'r', 'DisplayName', '$q_x$');
plot(T.timestamp_s, T.ATT_Q2, 'g', 'DisplayName', '$q_y$');
plot(T.timestamp_s, T.ATT_Q3, 'b', 'DisplayName', '$q_z$');
xlabel("Time [s]")
ylabel("Quaternion")
ylim([-1, 1])
grid on
% title("Attitude quaternion") 
legend('Location','southwest'); hold off;


f_e = figure(2);
plot(T.timestamp_s, T.euler_roll, 'r', 'DisplayName', '$\phi$'); hold on;
plot(T.timestamp_s, T.euler_pitch, 'g', 'DisplayName', '$\theta$');
plot(T.timestamp_s, T.euler_yaw, 'b', 'DisplayName', '$\psi$');
xlabel("Time [s]")
ylabel("Angle [rad]")
grid on
% title("Relative Euler angles")
legend('Location','southwest'); hold off;

f_w = figure(3);
plot(T.timestamp_s, T.RATE_WX, 'r', 'DisplayName', '$\omega_x$'); hold on;
plot(T.timestamp_s, T.RATE_WY, 'g', 'DisplayName', '$\omega_y$')
plot(T.timestamp_s, T.RATE_WZ, 'b', 'DisplayName', '$\omega_z$')
xlabel("Time [s]")
ylabel("Angular velocity [rad/s]")
grid on
% title("Angular rates")
legend('Location','northwest'); hold off;

f_v = figure(4);
plot(T.timestamp_s, T.VEL_VX, 'r', 'DisplayName', '$v_x$'); hold on;
plot(T.timestamp_s, T.VEL_VY, 'g', 'DisplayName', '$v_y$');
plot(T.timestamp_s, T.VEL_VZ, 'b', 'DisplayName', '$v_z$');
xlabel("Time [s]")
ylabel("Velocity [m/s]")
grid on
% title("Velocity")
legend('Location','best'); hold off;

f_a = figure(5);
plot(T.timestamp_s, T.ALT, 'b', 'DisplayName', '$r_x$')
xlabel("Time [s]")
ylabel("Altitude [m]")
grid on
% title("Altitude")
legend('Location','best'); 
hold off;

f_c = figure(6);
plot(T.timestamp_s, rad2deg(T.CANARD_ANGLE), 'r', 'DisplayName', '$\delta$'); hold on;
plot(T.timestamp_s, T.COEFF_CL, 'b', 'DisplayName', '$C_L$')
xlabel("Time [s]")
ylabel("Angle [deg], Coefficient")
ylim([-1,5])
grid on
% title("Canard")
legend('Location','best'); hold off;

set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')

fontsize(f_q, 12, "points")
fontsize(f_e, 12, "points")
fontsize(f_w, 12, "points")
fontsize(f_v, 12, "points")
fontsize(f_a, 12, "points")
fontsize(f_c, 12, "points")

%% export
exportgraphics(f_q, 'analysis/testflight/testflight_q.pdf')
exportgraphics(f_e, 'analysis/testflight/testflight_euler.pdf')
exportgraphics(f_w, 'analysis/testflight/testflight_w.pdf')
exportgraphics(f_v, 'analysis/testflight/testflight_v.pdf')
exportgraphics(f_a, 'analysis/testflight/testflight_alt.pdf')
exportgraphics(f_c, 'analysis/testflight/testflight_canard.pdf')