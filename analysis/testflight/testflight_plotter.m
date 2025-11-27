%% load
load("analysis\testflight\testflight_flight.mat")

%% plot

%%% Colors (same as plot_state)
colors(1,:) = [0, 0, 0];      % Black
colors(2,:) = [0.9, 0.3, 0.1];% Red
colors(3,:) = [0.2, 0.8, 0.1];% Green
colors(4,:) = [0.1, 0.2, 0.8];% Blue

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')

f_tf = figure(1);
tiledlayout(3,2,'TileSpacing','Compact');

% f_q = figure(1);
p_q = nexttile;
plot(T.timestamp_s, T.ATT_Q0,  'Color', colors(1,:), 'DisplayName', '$q_w$'); hold on;
plot(T.timestamp_s, T.ATT_Q1,  'Color', colors(2,:), 'DisplayName', '$q_x$');
plot(T.timestamp_s, T.ATT_Q2,  'Color', colors(3,:), 'DisplayName', '$q_y$');
plot(T.timestamp_s, T.ATT_Q3,  'Color', colors(4,:), 'DisplayName', '$q_z$');
xlabel("Time [s]")
% ylabel("Quaternion")
title("Quaternion")
ylim([-1, 1])
grid on
% title("Attitude quaternion") 
% legend('Location','southwest'); 
legend('Orientation','vertical', Location='eastoutside', IconColumnWidth=12)
hold off;


% f_e = figure(2);
nexttile;
plot(T.timestamp_s, T.euler_roll, 'Color', colors(2,:), 'DisplayName', '$\phi$'); hold on;
plot(T.timestamp_s, T.euler_pitch, 'Color', colors(3,:), 'DisplayName', '$\theta$');
plot(T.timestamp_s, T.euler_yaw, 'Color', colors(4,:), 'DisplayName', '$\psi$');
xlabel("Time [s]")
% ylabel("Angle [rad]")
title("Euler angle [rad]")
grid on
% title("Relative Euler angles")
% legend('Location','southwest'); 
legend('Orientation','vertical', Location='eastoutside', IconColumnWidth=12)
hold off;

% f_w = figure(3);
nexttile;
plot(T.timestamp_s, T.RATE_WX,  'Color', colors(2,:), 'DisplayName', '$\omega_x$'); hold on;
plot(T.timestamp_s, T.RATE_WY,  'Color', colors(3,:), 'DisplayName', '$\omega_y$')
plot(T.timestamp_s, T.RATE_WZ,  'Color', colors(4,:), 'DisplayName', '$\omega_z$')
xlabel("Time [s]")
% ylabel("Angular velocity [rad/s]")
title("Angular velocity [rad/s]")
grid on
% title("Angular rates")
% legend('Location','northwest'); 
legend('Orientation','vertical', Location='eastoutside', IconColumnWidth=12)
hold off;

% f_v = figure(4);
nexttile;
plot(T.timestamp_s, T.VEL_VX,  'Color', colors(2,:), 'DisplayName', '$v_x$'); hold on;
plot(T.timestamp_s, T.VEL_VY,  'Color', colors(3,:), 'DisplayName', '$v_y$');
plot(T.timestamp_s, T.VEL_VZ,  'Color', colors(4,:), 'DisplayName', '$v_z$');
xlabel("Time [s]")
% ylabel("Velocity [m/s]")
title("Velocity [m/s]")
grid on
% title("Velocity")
% legend('Location','best'); 
legend('Orientation','vertical', Location='eastoutside', IconColumnWidth=12)
hold off;

% f_a = figure(5);
nexttile;
plot(T.timestamp_s, T.ALT, 'b', 'DisplayName', '$r_x$')
xlabel("Time [s]")
% ylabel("Altitude [m]")
title("Altitude [m]")
grid on
% title("Altitude")
% legend('Location','best');
legend('Orientation','vertical', Location='eastoutside', IconColumnWidth=12)
hold off;

% f_c = figure(6);
nexttile;
plot(T.timestamp_s, rad2deg(T.CANARD_ANGLE),  'Color', colors(2,:), 'DisplayName', '$\delta$'); hold on;
plot(T.timestamp_s, T.COEFF_CL,  'Color', colors(4,:), 'DisplayName', '$C_L$')
xlabel("Time [s]")
% ylabel("Angle [deg], Coefficient")
title("Angle [deg], Coefficient")
ylim([-1,5])
grid on
% title("Canard")
% legend('Location','best'); 
legend('Orientation','vertical', Location='eastoutside', IconColumnWidth=12)
hold off;

f_descent = figure(2);
for t=1:height(T)
    q = table2array(T(t,2:5))';
    v = table2array(T(t,11:13))';
    S = quaternion_rotmatrix(q);
    v_G = S'*v;
    T.VEL_VERT(t) = v_G(1);
end
tiledlayout(2,1,'TileSpacing','Compact');
nexttile
plot(T.timestamp_s, smooth(T.VEL_VERT, 10));
hold on
plot(T.timestamp_s, smooth([0;diff(T.ALT)*15], 10));
hold off
nexttile
plot(T.timestamp_s, T.ALT) 
hold on
plot(T.timestamp_s, cumtrapz(smooth(T.VEL_VERT, 10)/15)+250);
hold off

set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')

fontsize(f_tf, 12, "points")
% fontsize(f_q, 12, "points")
% fontsize(f_e, 12, "points")
% fontsize(f_w, 12, "points")
% fontsize(f_v, 12, "points")
% fontsize(f_a, 12, "points")
% fontsize(f_c, 12, "points")

%% export
exportgraphics(f_tf, 'analysis/testflight/testflight.pdf')
% exportgraphics(f_q, 'analysis/testflight/testflight_q.pdf')
% exportgraphics(f_e, 'analysis/testflight/testflight_euler.pdf')
% exportgraphics(f_w, 'analysis/testflight/testflight_w.pdf')
% exportgraphics(f_v, 'analysis/testflight/testflight_v.pdf')
% exportgraphics(f_a, 'analysis/testflight/testflight_alt.pdf')
% exportgraphics(f_c, 'analysis/testflight/testflight_canard.pdf')