% Path to your data
datapath = "analysis/benchtests/test/";

% --- File groups ---
imu_files = [ ...
    "movella_pt1.csv", "movella_pt2.csv", "movella_pt3.csv", ...
    "pololu_pt1.csv",  "pololu_pt2.csv",  "pololu_pt3.csv" ];

ekf_control_files = [ ...
    "ekf_ctx_pt1.csv", "ekf_ctx_pt2.csv", "ekf_ctx_pt3.csv", ...
    "controller_input" ];

% --- Helper function: read table -> remove timestamp -> timetable ---
function tt = importCSV(fname)
    [~,base,~] = fileparts(fname);                 % filename without .csv
    tbl = readtable(fname,'PreserveVariableNames',true);

    % Convert ms → seconds, then round to nearest millisecond
    t = duration( round( seconds(tbl.timestamp/1000), 3 ), 0, 0 );  
    % now t is a duration vector, rounded to 0.001 s

    tbl.timestamp = [];                            % drop timestamp column

    % Prefix each variable with the filename
    tbl = renamevars(tbl, tbl.Properties.VariableNames, ...
                          strcat(base,"_",tbl.Properties.VariableNames));

    % Build timetable
    tt  = table2timetable(tbl,'RowTimes',t);
end


% --- Import and synchronize Movella + Pololu ---
TTs_MP = cellfun(@(f) importCSV(datapath+f), imu_files, 'UniformOutput', false);
imu_tt = synchronize(TTs_MP{:});

% --- Import and synchronize EKF + Control ---
TTs_EC = cellfun(@(f) importCSV(datapath+f), ekf_control_files, 'UniformOutput', false);
est_tt = synchronize(TTs_EC{:});

%%

T_pad_start = 24820 / 1000;
T_efk_start = 49195 / 1000;

sensor_select = [0,1,1]';

xhat = zeros(13,1); Phat = zeros(13); bias_1 = zeros(10, 1); bias_2 = zeros(10, 1);
x = xhat; P = Phat; b.bias_1 = bias_1; b.bias_2 = bias_2;
x_array = nan(13, length(timestamps));
x_array(:,1) = xhat;

t = timestamps(1);
for k = 2:length(timestamps)

    dt = timestamps(k) - t;
    t = timestamps(k);

    IMU_2 = IMU_array(k,:)';
    
    if t >= T_pad_start && t < T_efk_start
        [xhat, bias_1, bias_2] = pad_filter(IMU_1, IMU_2, sensor_select(1:2));
        x = xhat; b.bias_1 = bias_1; b.bias_2 = bias_2;
    end
    if t >= T_efk_start
        [xhat, Phat] = ekf_algorithm(x, P, b, t, dt, IMU_1, IMU_2, cmd, encoder, sensor_select);
        x = xhat; P = Phat;
    end
    x_array(:,k) = x;
end

x_array(:,end)

% %% Plot

% figure(1)
% plot(timestamps, x_array(1:4,:)');
% xlabel("Time [s]")
% legend(["w", "x", "y", "z"]);

% %%
% 
% figure(2)
% stairs(timestamps, T.euler_roll, 'DisplayName', 'roll'); hold on;
% stairs(timestamps, T.euler_pitch, 'DisplayName', 'pitch');
% stairs(timestamps, T.euler_yaw, 'DisplayName', 'yaw');
% xlabel("Time [s]")
% ylabel("Angle [rad]")
% % title("Relative Euler angles")
% legend('Location','southwest'); hold off;
% 
% figure(3)
% stairs(timestamps, T.RATE_WX, 'DisplayName', '$\omega_x$'); hold on;
% stairs(timestamps, T.RATE_WY, 'DisplayName', '$\omega_y$')
% stairs(timestamps, T.RATE_WZ, 'DisplayName', '$\omega_z$')
% xlabel("Time [s]")
% ylabel("Anglular rate [rad/s]")
% % title("Angular rates")
% legend('Location','northwest'); hold off;
% 
% figure(4)
% stairs(timestamps, T.VEL_VX, 'DisplayName', '$v_x$'); hold on;
% stairs(timestamps, T.VEL_VY, 'DisplayName', '$v_y$');
% stairs(timestamps, T.VEL_VZ, 'DisplayName', '$v_z$');
% xlabel("Time [s]")
% ylabel("Velocity [m/s]")
% % title("Velocity")
% legend('Location','best'); hold off;
% 
% figure(5)
% stairs(timestamps, T.ALT, 'DisplayName', 'alt')
% xlabel("Time [s]")
% ylabel("Altitude [m]")
% % title("Altitude")
% %legend(); 
% hold off;
% 
% figure(6)
% stairs(timestamps, rad2deg(T.CANARD_ANGLE), 'DisplayName', '$\delta$'); hold on;
% stairs(timestamps, T.COEFF_CL, 'DisplayName', '$C_L$')
% xlabel("Time [s]")
% ylabel("Angle [deg], Coefficient [ ]")
% ylim([-1,5])
% % title("Canard")
% legend('Location','best'); hold off;