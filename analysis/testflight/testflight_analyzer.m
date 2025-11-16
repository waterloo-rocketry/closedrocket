clear
% Read data (if not already in workspace)
T_long = readtable('analysis/testflight/testflightomnibusloganal.xlsx', 'Sheet', 'proc_state_est'); 

% relevant times only
% keep only times between time_start and time_end
time_start = 1455; % ekf starts at 1403, liftoff around 1463
time_end = 1486;
time_end = 1500;
time_end = time_start + 40;

T_long.Var4 = [];
T_long.notes = [];

% Convert long to wide format
T = unstack(T_long, 'data', 'state_id', AggregationFunction=@mean);
% format names
oldNames = T.Properties.VariableNames;
newNames = strrep(oldNames, 'STATE_ID_', '');
T.Properties.VariableNames = newNames;

% replace NaN with previous 
T = fillmissing(T, 'linear');

% % Convert from milliseconds to seconds
T.timestamp_s = T.timestamp_ms / 1000;
T.timestamp_ms = [];

% rebase to new zero time
T = T(T.timestamp_s >= time_start & T.timestamp_s <= time_end, :);
T.timestamp_s = T.timestamp_s - T.timestamp_s(1);


%% process data
% euler angles
for i=1:height(T)
    q = [T.ATT_Q0(i), T.ATT_Q1(i), T.ATT_Q2(i), T.ATT_Q3(i)]';
    euler = quaternion_to_euler(q);
    T.euler_roll(i) = euler(1);
    T.euler_pitch(i) = euler(2);
    T.euler_yaw(i) = euler(3);
end

save("analysis\testflight\testflight_flight.mat", "T")

