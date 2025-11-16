clear
%% Read files (if not already in workspace)
T1_est = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'state_est_data'); 
T1_cmd = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'proc_cmd'); 
T1_enc = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'mcb_encoder'); 
T1_imu = readtable('analysis/aurora/aurora_flight_data.xlsx', 'Sheet', 'altimu_meas'); 

T2_est = readtable('analysis/aurora/aurora_estimator_controller_complete.xlsx', 'Sheet', 'state_est_data'); 
T2_enc = readtable('analysis/aurora/aurora_estimator_controller_complete.xlsx', 'Sheet', 'mcb_encoder'); 
T2_imu = readtable('analysis/aurora/aurora_estimator_controller_complete.xlsx', 'Sheet', 'altimu_meas'); 

%% unpack tables
% Convert long to wide format
T1 = unstack(T1_est, 'data', 'state_id', AggregationFunction=@mean);
T2 = unstack(T2_est, 'data', 'state_id', AggregationFunction=@mean);
% format names
oldNames = T1.Properties.VariableNames;
newNames = strrep(oldNames, 'STATE_ID_', '');
T1.Properties.VariableNames = newNames;
oldNames = T2.Properties.VariableNames;
newNames = strrep(oldNames, 'STATE_ID_', '');
T2.Properties.VariableNames = newNames;


%% synchronize
% relevant times only
% keep only times between time_start and time_end
time1_start = 0;
time1_end = 1500;
time1_proc_offset = 10863;
time1_mcb_offset = 30644;% -2e3;

time2_start = 0;
time2_end = 1500;
time2_proc_offset = 169219983 + 592000;
time2_mcb_offset = 656625343;% -2e3;

% Convert from milliseconds to seconds
function T = retimer(T, time_offset, time_start, time_end)
    
    T.time_s = seconds((T.time_ms - time_offset)/ 1000 );
    T.time_ms = [];
    T = T(T.time_s >= time_start & T.time_s <= time_end, :);
    % T.time_s = T.time_s - T.time_s(1);
end

T1 = table2timetable(retimer(T1, time1_proc_offset, time1_start, time1_end), 'RowTimes','time_s');
T1_imu = table2timetable(retimer(T1_imu, time1_proc_offset, time1_start, time1_end), 'RowTimes','time_s');
T1_cmd = table2timetable(retimer(T1_cmd, time1_proc_offset, time1_start, time1_end), 'RowTimes','time_s');
T1_enc = table2timetable(retimer(T1_enc, time1_mcb_offset, time1_start, time1_end), 'RowTimes','time_s');

T2 = table2timetable(retimer(T2, time2_proc_offset, time2_start, time2_end), 'RowTimes','time_s');
T2_imu = table2timetable(retimer(T2_imu, time2_proc_offset, time2_start, time2_end), 'RowTimes','time_s');
T2_enc = table2timetable(retimer(T2_enc, time2_mcb_offset, time2_start, time2_end), 'RowTimes','time_s');


T = sortrows([T1; T2]);
T_enc = sortrows([T1_enc; T2_enc]);
% T_imu = sortrows([T1_imu; T2_imu]);
T_imu = T1_imu; 
T_cmd = T1_cmd;


%% process data

% replace NaN with interpolation
TF = fillmissing(T, 'previous');
TF_imu = fillmissing(T_imu, 'previous');
% TF_enc = fillmissing(T_enc, 'previous');

% TL_enc = fillmissing(T_enc, 'linear');



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


%%% Interpolation
TL = fillmissing(T, 'linear', 'EndValues', 'none');
TL_imu = fillmissing(T_imu, 'linear', 'EndValues', 'none');


%%% euler angles
for i=1:height(T)
    q = [TL.ATT_Q0(i), TL.ATT_Q1(i), TL.ATT_Q2(i), TL.ATT_Q3(i)]';
    euler = quaternion_to_euler(q);
    TL.euler_roll(i) = euler(1);
    TL.euler_pitch(i) = euler(2);
    TL.euler_yaw(i) = euler(3);
end

%% print Recovery
fprintf("Recovery: at %f s and %f m altitude \n the velocity was %f m/s vertical, %f m/s lateral.\n", seconds(T.time_s(12)), T.ALT(12), vel_apo_vert, vel_apo_hor);


%% save
save("analysis\aurora\aurora_flight.mat", "T", "TL", "T_imu", "TL_imu", "T_cmd", "T_enc");

