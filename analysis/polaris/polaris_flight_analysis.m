%% Polaris flight data analysis
% Parse the sparse flight-computer CSV and display navigation, controller,
% sensor, and actuator dashboards. Each telemetry message is kept at its
% native sample times; missing cells are not interpolated or forward-filled.

%% Configuration
script_dir = fileparts(mfilename("fullpath"));
csv_file = fullfile(script_dir, "00000048.BIN-parsed.csv");

% Time-zero reference. The first controller record is about 7 seconds after
% liftoff in this log, so subtract that delay to put liftoff at t = 0.
% time_zero_source may also be "navigation", "log", or an absolute CSV
% timestamp in milliseconds.
time_zero_source = "controller";
controller_start_after_liftoff = 7; % s

% Plot window in seconds after liftoff. Set both to [] for the complete log.
plot_start_time = 0;
plot_end_time = 30;
if isempty(plot_start_time) && isempty(plot_end_time)
    time_limits = [];
else
    time_limits = [plot_start_time, plot_end_time];
end

save_figures = true;
output_dir = fullfile(script_dir, "figures");

%% Parse and plot
polaris = parse_polaris_csv(csv_file, time_zero_source, ...
    controller_start_after_liftoff);
figs = plot_polaris_data(polaris, time_limits);
%%
if save_figures
    save_polaris_figures(figs, output_dir);
end

%% Local functions
function data = parse_polaris_csv(csvFile, timeZero, controllerDelay)
    if ~isfile(csvFile)
        error("polaris:FileNotFound", "CSV file not found: %s", csvFile);
    end

    fprintf("Reading Polaris flight data from %s ...\n", csvFile);
    raw = readtable(csvFile, "VariableNamingRule", "preserve");

    required = [ ...
        "timestamp", ...
        "navigator_pt1.q_w", "navigator_pt1.q_x", ...
        "navigator_pt1.q_y", "navigator_pt1.q_z", ...
        "navigator_pt1.altitude", "navigator_pt1.norm", ...
        "navigator_pt2.vel_x", "navigator_pt2.vel_y", ...
        "navigator_pt2.vel_z", "navigator_pt2.rate_x", ...
        "navigator_pt2.rate_y", "navigator_pt2.rate_z", ...
        "controller.command", "controller.canard_coeff", ...
        "controller.body_coeff", "controller.roll_angle_target", ...
        "controller.roll_rate_target", ...
        "board_imu.accel_x", "board_imu.accel_y", "board_imu.accel_z", ...
        "board_imu.gyro_x", "board_imu.gyro_y", "board_imu.gyro_z", ...
        "board_mag_baro.mag_x", "board_mag_baro.mag_y", ...
        "board_mag_baro.mag_z", "board_mag_baro.baro", ...
        "board_mag_baro.temp", ...
        "movella_pt1.accel_x", "movella_pt1.accel_y", ...
        "movella_pt1.accel_z", "movella_pt1.gyro_x", ...
        "movella_pt1.gyro_y", "movella_pt1.gyro_z", ...
        "movella_pt2.mag_x", "movella_pt2.mag_y", ...
        "movella_pt2.mag_z", "movella_pt2.baro", ...
        "movella_pt3.q_w", "movella_pt3.q_x", ...
        "movella_pt3.q_y", "movella_pt3.q_z", ...
        "ad_breakout.accel_x", "ad_breakout.accel_y", ...
        "ad_breakout.accel_z", "ad_breakout.gyro", ...
        "servo_motor.angle", "servo_motor.curr", "servo_motor.temp"];

    available = string(raw.Properties.VariableNames);
    missing = required(~ismember(required, available));
    if ~isempty(missing)
        error("polaris:MissingColumns", ...
            "CSV is missing required columns: %s", strjoin(missing, ", "));
    end

    timeZeroMs = resolve_time_zero(raw, timeZero);
    if isstring(timeZero) || ischar(timeZero)
        if lower(string(timeZero)) == "controller"
            timeZeroMs = timeZeroMs - 1000 * controllerDelay;
        end
    end

    navPt1Columns = [ ...
        "navigator_pt1.q_w", "navigator_pt1.q_x", ...
        "navigator_pt1.q_y", "navigator_pt1.q_z", ...
        "navigator_pt1.altitude", "navigator_pt1.norm"];
    rows = populated_rows(raw, navPt1Columns);
    data.navigation.time_pt1 = relative_time(raw, rows, timeZeroMs);
    data.navigation.q = values(raw, rows, navPt1Columns(1:4));
    data.navigation.altitude = values(raw, rows, navPt1Columns(5));
    data.navigation.covariance_norm = values(raw, rows, navPt1Columns(6));
    data.navigation.roll = quaternion_roll(data.navigation.q);

    navPt2Columns = [ ...
        "navigator_pt2.vel_x", "navigator_pt2.vel_y", ...
        "navigator_pt2.vel_z", "navigator_pt2.rate_x", ...
        "navigator_pt2.rate_y", "navigator_pt2.rate_z"];
    rows = populated_rows(raw, navPt2Columns);
    data.navigation.time_pt2 = relative_time(raw, rows, timeZeroMs);
    data.navigation.velocity = values(raw, rows, navPt2Columns(1:3));
    data.navigation.rates = values(raw, rows, navPt2Columns(4:6));

    controllerColumns = [ ...
        "controller.command", "controller.canard_coeff", ...
        "controller.body_coeff", "controller.roll_angle_target", ...
        "controller.roll_rate_target"];
    rows = populated_rows(raw, controllerColumns);
    data.controller.time = relative_time(raw, rows, timeZeroMs);
    data.controller.command = values(raw, rows, controllerColumns(1));
    data.controller.canard_coefficient = values(raw, rows, controllerColumns(2));
    data.controller.body_coefficient = values(raw, rows, controllerColumns(3));
    data.controller.roll_angle_target = values(raw, rows, controllerColumns(4));
    data.controller.roll_rate_target = values(raw, rows, controllerColumns(5));

    boardImuColumns = [ ...
        "board_imu.accel_x", "board_imu.accel_y", "board_imu.accel_z", ...
        "board_imu.gyro_x", "board_imu.gyro_y", "board_imu.gyro_z"];
    rows = populated_rows(raw, boardImuColumns);
    data.board_imu.time = relative_time(raw, rows, timeZeroMs);
    data.board_imu.acceleration = values(raw, rows, boardImuColumns(1:3));
    data.board_imu.rates = values(raw, rows, boardImuColumns(4:6));

    boardEnvironmentColumns = [ ...
        "board_mag_baro.mag_x", "board_mag_baro.mag_y", ...
        "board_mag_baro.mag_z", "board_mag_baro.baro", ...
        "board_mag_baro.temp"];
    rows = populated_rows(raw, boardEnvironmentColumns);
    data.board_environment.time = relative_time(raw, rows, timeZeroMs);
    data.board_environment.magnetic_field = ...
        values(raw, rows, boardEnvironmentColumns(1:3));
    data.board_environment.pressure = values(raw, rows, boardEnvironmentColumns(4));
    data.board_environment.temperature = values(raw, rows, boardEnvironmentColumns(5));

    movellaImuColumns = [ ...
        "movella_pt1.accel_x", "movella_pt1.accel_y", ...
        "movella_pt1.accel_z", "movella_pt1.gyro_x", ...
        "movella_pt1.gyro_y", "movella_pt1.gyro_z"];
    rows = populated_rows(raw, movellaImuColumns);
    data.movella_imu.time = relative_time(raw, rows, timeZeroMs);
    data.movella_imu.acceleration = values(raw, rows, movellaImuColumns(1:3));
    data.movella_imu.rates = values(raw, rows, movellaImuColumns(4:6));

    movellaEnvironmentColumns = [ ...
        "movella_pt2.mag_x", "movella_pt2.mag_y", ...
        "movella_pt2.mag_z", "movella_pt2.baro"];
    rows = populated_rows(raw, movellaEnvironmentColumns);
    data.movella_environment.time = relative_time(raw, rows, timeZeroMs);
    data.movella_environment.magnetic_field = ...
        values(raw, rows, movellaEnvironmentColumns(1:3));
    data.movella_environment.pressure = ...
        values(raw, rows, movellaEnvironmentColumns(4));

    movellaAttitudeColumns = [ ...
        "movella_pt3.q_w", "movella_pt3.q_x", ...
        "movella_pt3.q_y", "movella_pt3.q_z"];
    rows = populated_rows(raw, movellaAttitudeColumns);
    data.movella_attitude.time = relative_time(raw, rows, timeZeroMs);
    data.movella_attitude.q = values(raw, rows, movellaAttitudeColumns);

    adColumns = [ ...
        "ad_breakout.accel_x", "ad_breakout.accel_y", ...
        "ad_breakout.accel_z", "ad_breakout.gyro"];
    rows = populated_rows(raw, adColumns);
    data.ad_breakout.time = relative_time(raw, rows, timeZeroMs);
    data.ad_breakout.acceleration = values(raw, rows, adColumns(1:3));
    data.ad_breakout.rate = values(raw, rows, adColumns(4));

    servoColumns = [ ...
        "servo_motor.angle", "servo_motor.curr", "servo_motor.temp"];
    rows = populated_rows(raw, servoColumns);
    data.servo.time = relative_time(raw, rows, timeZeroMs);
    data.servo.angle = values(raw, rows, servoColumns(1));
    data.servo.current = values(raw, rows, servoColumns(2));
    data.servo.temperature = values(raw, rows, servoColumns(3));

    data.meta.csv_file = string(csvFile);
    data.meta.time_zero_ms = timeZeroMs;
    data.meta.time_zero_source = string(timeZero);
    data.meta.controller_start_after_liftoff = controllerDelay;
    fprintf("Parsed %d CSV rows; plot time zero is %.0f ms.\n", ...
        height(raw), timeZeroMs);
end

function timeZeroMs = resolve_time_zero(raw, source)
    if isnumeric(source) && isscalar(source) && isfinite(source)
        timeZeroMs = double(source);
        return;
    end

    source = lower(string(source));
    switch source
        case "controller"
            column = "controller.command";
        case "navigation"
            column = "navigator_pt1.q_w";
        case "log"
            timeZeroMs = raw.timestamp(1);
            return;
        otherwise
            error("polaris:InvalidTimeZero", ...
                "TimeZero must be 'controller', 'navigation', 'log', or a numeric timestamp in ms.");
    end

    first = find(isfinite(raw.(column)), 1, "first");
    if isempty(first)
        error("polaris:EmptyTimeSource", ...
            "No samples are available for the '%s' time-zero source.", source);
    end
    timeZeroMs = raw.timestamp(first);
end

function rows = populated_rows(raw, columns)
    rows = any(isfinite(raw{:, cellstr(columns)}), 2);
end

function result = values(raw, rows, columns)
    result = raw{rows, cellstr(columns)};
end

function time = relative_time(raw, rows, timeZeroMs)
    time = (raw.timestamp(rows) - timeZeroMs) / 1000;
end

function roll = quaternion_roll(q)
    qw = q(:, 1);
    qx = q(:, 2);
    qy = q(:, 3);
    qz = q(:, 4);
    roll = atan2(2 .* (qy .* qz + qw .* qx), ...
        qw.^2 - qx.^2 - qy.^2 + qz.^2);
end

function figs = plot_polaris_data(data, timeLimits)
    colors = component_colors();
    figs = struct();

    figs.navigation = tabbed_figure("Navigation State Estimates");
    layout = tiledlayout(figs.navigation, 2, 2, ...
        "TileSpacing", "compact", "Padding", "compact");
    ax = nexttile(layout);
    plot_matrix(ax, data.navigation.time_pt1, data.navigation.q, ...
        ["q_w", "q_x", "q_y", "q_z"], colors.wxyz, "-");
    finish_axis(ax, "Attitude Quaternion", "", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.navigation.time_pt2, data.navigation.velocity, ...
        ["v_x", "v_y", "v_z"], colors.xyz, "-");
    finish_axis(ax, "Body-Frame Velocity", "m/s", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.navigation.time_pt2, data.navigation.rates, ...
        ["omega_x", "omega_y", "omega_z"], colors.xyz, "-");
    finish_axis(ax, "Angular Rates", "rad/s", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.navigation.time_pt1, data.navigation.altitude, ...
        "altitude", colors.z, "-");
    finish_axis(ax, "Altitude", "m", timeLimits);

    figs.covariance = tabbed_figure("Navigation Covariance");
    ax = axes(figs.covariance);
    plot_matrix(ax, data.navigation.time_pt1, data.navigation.covariance_norm, ...
        "covariance_norm", colors.z, "-");
    finish_axis(ax, "Covariance Norm", "", timeLimits);

    figs.controller = tabbed_figure("Controller Dashboard");
    layout = tiledlayout(figs.controller, 2, 2, ...
        "TileSpacing", "compact", "Padding", "compact");
    ax = nexttile(layout);
    plot_matrix(ax, data.controller.time, data.controller.roll_angle_target, ...
        "phi_target", colors.x, "-");
    plot_matrix(ax, data.navigation.time_pt1, data.navigation.roll, ...
        "phi", colors.z, "-");
    finish_axis(ax, "Roll Angle", "rad", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.controller.time, rad2deg(data.controller.command), ...
        "command", colors.x, "-");
    plot_matrix(ax, data.servo.time, data.servo.angle, ...
        "servo_angle", colors.z, "-");
    finish_axis(ax, "Motor Angle", "deg", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.controller.time, data.controller.roll_rate_target, ...
        "omega_target", colors.x, "-");
    plot_matrix(ax, data.navigation.time_pt2, data.navigation.rates(:, 1), ...
        "omega_x", colors.z, "-");
    finish_axis(ax, "Roll Rate", "rad/s", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.controller.time, ...
        [data.controller.canard_coefficient, data.controller.body_coefficient], ...
        ["C_l_delta", "C_l_0"], [colors.x; colors.z], "-");
    finish_axis(ax, "Estimated Roll Coefficients", "", timeLimits);

    figs.imus = tabbed_figure("Sensor IMUs");
    layout = tiledlayout(figs.imus, 2, 3, ...
        "TileSpacing", "compact", "Padding", "compact");
    ax = nexttile(layout);
    plot_matrix(ax, data.board_imu.time, data.board_imu.acceleration, ...
        ["a_x", "a_y", "a_z"], colors.xyz, "-");
    finish_axis(ax, "Board Accelerometer", "m/s^2", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.movella_imu.time, data.movella_imu.acceleration, ...
        ["a_x", "a_y", "a_z"], colors.xyz, "-");
    finish_axis(ax, "Movella Accelerometer", "m/s^2", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.ad_breakout.time, data.ad_breakout.acceleration, ...
        ["a_x", "a_y", "a_z"], colors.xyz, "-");
    finish_axis(ax, "AD Breakout Accelerometer", "m/s^2", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.board_imu.time, data.board_imu.rates, ...
        ["omega_x", "omega_y", "omega_z"], colors.xyz, "-");
    finish_axis(ax, "Board Gyroscope", "rad/s", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.movella_imu.time, data.movella_imu.rates, ...
        ["omega_x", "omega_y", "omega_z"], colors.xyz, "-");
    finish_axis(ax, "Movella Gyroscope", "rad/s", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.ad_breakout.time, data.ad_breakout.rate, ...
        "omega", colors.x, "-");
    finish_axis(ax, "AD Breakout Gyroscope", "rad/s", timeLimits);

    figs.environment = tabbed_figure("Environmental Sensors");
    layout = tiledlayout(figs.environment, 2, 2, ...
        "TileSpacing", "compact", "Padding", "compact");
    ax = nexttile(layout);
    plot_matrix(ax, data.board_environment.time, ...
        data.board_environment.magnetic_field, ...
        ["mag_x", "mag_y", "mag_z"], colors.xyz, "-");
    finish_axis(ax, "Board Magnetometer", "", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.movella_environment.time, ...
        data.movella_environment.magnetic_field, ...
        ["mag_x", "mag_y", "mag_z"], colors.xyz, "-");
    finish_axis(ax, "Movella Magnetometer", "", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.board_environment.time, ...
        data.board_environment.pressure / 1000, ...
        "board_pressure", colors.x, "-");
    plot_matrix(ax, data.movella_environment.time, ...
        data.movella_environment.pressure / 1000, ...
        "movella_pressure", colors.z, "-");
    finish_axis(ax, "Barometric Pressure", "kPa", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.movella_attitude.time, data.movella_attitude.q, ...
        ["q_w", "q_x", "q_y", "q_z"], colors.wxyz, "-");
    finish_axis(ax, "Movella Attitude Quaternion", "", timeLimits);

    figs.servo = tabbed_figure("Servo Motor Telemetry");
    layout = tiledlayout(figs.servo, 3, 1, ...
        "TileSpacing", "compact", "Padding", "compact");
    ax = nexttile(layout);
    plot_matrix(ax, data.servo.time, data.servo.angle, ...
        "servo_angle", colors.z, "-");
    finish_axis(ax, "Motor Angle", "deg", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.servo.time, data.servo.current, ...
        "servo_current", colors.x, "-");
    finish_axis(ax, "Motor Current", "A", timeLimits);
    ax = nexttile(layout);
    plot_matrix(ax, data.servo.time, data.servo.temperature, ...
        "servo_temperature", colors.y, "-");
    finish_axis(ax, "Motor Temperature", "deg C", timeLimits);
end

function plot_matrix(ax, time, matrix, labels, colors, lineStyle)
    hold(ax, "on");
    labels = string(labels);
    if isvector(matrix)
        matrix = matrix(:);
    end
    if size(colors, 1) == 1 && size(matrix, 2) > 1
        colors = repmat(colors, size(matrix, 2), 1);
    end
    for column = 1:size(matrix, 2)
        line = plot(ax, time, matrix(:, column), lineStyle, ...
            "LineWidth", 1.2, "DisplayName", labels(column));
        if size(colors, 1) >= column
            line.Color = colors(column, :);
        end
    end
end

function finish_axis(ax, titleText, yLabel, timeLimits)
    title(ax, titleText, "FontWeight", "normal", "Interpreter", "none");
    xlabel(ax, "Time (s)");
    ylabel(ax, yLabel);
    grid(ax, "on");
    box(ax, "on");
    ax.XMinorTick = "on";
    ax.YMinorTick = "on";
    ax.LineWidth = 0.8;
    if ~isempty(timeLimits)
        xlim(ax, timeLimits(1:2));
    end
    plotLegend = legend(ax, "show", "Location", "best", ...
        "NumColumns", 2, "Interpreter", "none");
    plotLegend.ItemTokenSize = [12, 18];
end

function colors = component_colors()
    colors.x = [0.8500, 0.3250, 0.0980];
    colors.y = [0.4660, 0.6740, 0.1880];
    colors.z = [0.0000, 0.4470, 0.7410];
    colors.w = [0.1000, 0.1000, 0.1000];
    colors.xyz = [colors.x; colors.y; colors.z];
    colors.wxyz = [colors.w; colors.xyz];
end

function fig = tabbed_figure(name)
    if usejava("desktop")
        fig = figure("Name", name, "WindowStyle", "docked");
    else
        % MATLAB batch sessions have no desktop in which to dock figures.
        fig = figure("Name", name);
    end
end

function save_polaris_figures(figs, outputDir)
    if ~isfolder(outputDir)
        mkdir(outputDir);
    end
    names = fieldnames(figs);
    screenSize = get(groot, "ScreenSize");
    screenDpi = get(groot, "ScreenPixelsPerInch");
    exportWidth = screenSize(3) / screenDpi;
    exportHeight = screenSize(4) / screenDpi;
    for index = 1:numel(names)
        drawnow;
        exportgraphics(figs.(names{index}), ...
            fullfile(outputDir, string(names{index}) + ".pdf"), ...
            "ContentType", "vector", ...
            "Width", exportWidth, "Height", exportHeight, ...
            "Units", "inches", "Padding", "figure", ...
            "PreserveAspectRatio", "off");
    end
    fprintf("Saved %d figures to %s.\n", numel(names), outputDir);
end
