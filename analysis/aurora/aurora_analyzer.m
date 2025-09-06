% Read data (if not already in workspace)
% T_long = readtable('analysis/aurora/aurora_flight_data.xlsx'); 
T_est = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'state_est_data'); 
T_cmd = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'proc_cmd'); 
T_enc = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'mcb_encoder'); 
%%

% relevant times only
% keep only times between time_start and time_end
time_start = 0;
time_end = 150;
time_proc_offset = 10863;
time_mcb_offset = 30644 -2e3;

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


%% plot est

f_q = figure(1);
plot(T.time_s, T.ATT_Q0, 'o:', 'DisplayName', 'w'); hold on;
plot(T.time_s, T.ATT_Q1, 'o:', 'DisplayName', 'x');
plot(T.time_s, T.ATT_Q2, 'o:', 'DisplayName', 'y');
plot(T.time_s, T.ATT_Q3, 'o:', 'DisplayName', 'z');
% plot(TF.time_s, TF.ATT_Q0, 'o:', 'DisplayName', 'w'); hold on;
% plot(TF.time_s, TF.ATT_Q1, 'o:', 'DisplayName', 'x');
% plot(TF.time_s, TF.ATT_Q2, 'o:', 'DisplayName', 'y');
% plot(TF.time_s, TF.ATT_Q3, 'o:', 'DisplayName', 'z');
xlabel("Time [s]")
ylabel("Quaternion [ ]")
ylim([-1, 1])
% title("Attitude quaternion") 
legend('Location','southeast'); hold off;

f_e = figure(2);
plot(T.time_s, T.euler_roll, 'o:', 'DisplayName', 'roll'); hold on;
plot(T.time_s, T.euler_pitch, 'o:', 'DisplayName', 'pitch');
plot(T.time_s, T.euler_yaw, 'o:', 'DisplayName', 'yaw');
xlabel("Time [s]")
ylabel("Angle [rad]")
% title("Relative Euler angles")
legend('Location','southwest'); hold off;

f_w = figure(3);
plot(T.time_s, T.RATE_WX, 'o:', 'DisplayName', 'x'); hold on;
plot(T.time_s, T.RATE_WY, 'o:', 'DisplayName', 'y')
plot(T.time_s, T.RATE_WZ, 'o:', 'DisplayName', 'z')
xlabel("Time [s]")
ylabel("Anglular rate [rad/s]")
% title("Angular rates")
legend('Location','northwest'); hold off;

f_v = figure(4);
plot(T.time_s, T.VEL_VX, 'o:', 'DisplayName', 'x'); hold on;
plot(T.time_s, T.VEL_VY, 'o:', 'DisplayName', 'y');
plot(T.time_s, T.VEL_VZ, 'o:', 'DisplayName', 'z');
xlabel("Time [s]")
ylabel("Velocity [m/s]")
% title("Velocity")
legend('Location','best'); hold off;

f_a = figure(5);
plot(T.time_s, T.ALT, 'o:', 'DisplayName', 'alt')
xlabel("Time [s]")
ylabel("Altitude [m]")
% title("Altitude")
%legend(); 
hold off;

f_c = figure(6);
plot(T.time_s, rad2deg(T.CANARD_ANGLE), 'o:', 'DisplayName', '\delta'); hold on;
plot(T.time_s, T.COEFF_CL, 'o:', 'DisplayName', 'C_L')
xlabel("Time [s]")
ylabel("Angle [deg], Coefficient [ ]")
% ylim([-1,5])
% title("Canard")
legend('Location','best'); hold off;

%% plot control

f_cmd = figure(7);
plot(T_cmd.time_s, T_cmd.data, '.:', 'DisplayName', 'cmd'); hold on;
plot(T_enc.time_s, T_enc.data, 'o-', 'DisplayName', 'enc')
xlabel("Time [s]")
ylabel("Command [deg], Encoder [deg]")
ylim([-12,12])
% title("Canard")
legend('Location','best'); hold off;


%% export
% exportgraphics(f_q, 'analysis/testflight/testflight_q.png')
% exportgraphics(f_e, 'analysis/testflight/testflight_euler.png')
% exportgraphics(f_w, 'analysis/testflight/testflight_w.png')
% exportgraphics(f_v, 'analysis/testflight/testflight_v.png')
% exportgraphics(f_a, 'analysis/testflight/testflight_alt.png')
% exportgraphics(f_c, 'analysis/testflight/testflight_canard.png')