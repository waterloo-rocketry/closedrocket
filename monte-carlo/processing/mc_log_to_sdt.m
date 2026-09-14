function sdt = mc_log_to_sdt(log)
%MC_LOG_TO_SDT Convert a closedrocket log to Monte Carlo plot timetables.

    validate_log(log);

    ref = signal_timetable(log, ["ref", "phi_target"], "ref");
    cmd = signal_timetable(log, "cmd", "cmd");
    roll = signal_timetable(log, ["roll", "phi"], "roll");
    delta = signal_timetable(log, ...
        ["delta_ctrl", "canard_angle", "delta"], "delta");
    w = signal_timetable(log, "w", "w");

    rate = [];
    if ~isempty(w)
        w_data = table2array(w);
        rate = timetable(w.Time, w_data(:, 1), ...
            'VariableNames', {'rate'});
    end

    roll_error = timetable_difference(ref, roll, "err");
    sdt.control = combine_timetables( ...
        {ref, roll_error, cmd, roll, rate, delta});

    q = signal_timetable(log, "q", "q");
    v = signal_timetable(log, "v", "v");
    pos = signal_timetable(log, ["pos", "r"], "pos");
    pos_yz = [];
    if ~isempty(pos)
        pos_data = table2array(pos);
        alt = timetable(pos.Time, pos_data(:, 1), ...
            'VariableNames', {'alt'});
        if size(pos_data, 2) >= 3
            pos_yz = timetable(pos.Time, pos_data(:, 2:3), ...
                'VariableNames', {'pos_yz'});
        end
    else
        alt = signal_timetable(log, ...
            ["alt", "alt_est", "alt_hat"], "alt");
    end

    cl = signal_timetable(log, ["cl", "CL", "CL_est"], "cl");
    rocket = combine_timetables({q, w, v, alt, cl, delta, pos_yz});

    q_hat = signal_timetable(log, ["q_hat", "q_est"], "q_hat");
    w_hat = signal_timetable(log, ["w_hat", "w_est"], "w_hat");
    v_hat = signal_timetable(log, ["v_hat", "v_est"], "v_hat");
    alt_hat = signal_timetable(log, ...
        ["alt_hat", "alt_est"], "alt_hat");
    sdt.est = combine_timetables({q_hat, w_hat, v_hat, alt_hat});

    if istimetable(sdt.est) && istimetable(rocket)
        sdt.rocket_dt = retime(rocket, sdt.est.Time, 'linear');
    else
        sdt.rocket_dt = rocket;
    end

    q_error = timetable_difference(q, q_hat, "q");
    w_error = timetable_difference(w, w_hat, "w");
    v_error = timetable_difference(v, v_hat, "v");
    alt_error = timetable_difference(alt, alt_hat, "alt");
    sdt.error = combine_timetables( ...
        {q_error, w_error, v_error, alt_error});

    p_norm = signal_timetable(log, "P_norm", "P_norm");
    if ~isempty(p_norm)
        sdt.P_norm = p_norm;
    end
end

function validate_log(log)
    if ~isstruct(log) || ~isscalar(log) || ...
            ~isfield(log, "format") || ...
            string(log.format) ~= "closedrocket-log" || ...
            ~isfield(log, "time") || ...
            ~isfield(log, "signals") || ...
            ~isstruct(log.signals)
        error("mc_log_to_sdt:InvalidLog", ...
            "Input must be a scalar closedrocket-log struct.");
    end
end

function tt = signal_timetable(log, aliases, variable_name)
    aliases = string(aliases);
    tt = [];
    for i = 1:numel(aliases)
        field = matlab.lang.makeValidName(char(aliases(i)));
        if ~isfield(log.signals, field)
            continue;
        end

        time = log.time(:);
        values = log.signals.(field);
        values = reshape(values, numel(time), []);
        tt = timetable(seconds(time), values, ...
            'VariableNames', {char(variable_name)});
        return;
    end
end

function combined = combine_timetables(tables)
    tables = tables(~cellfun(@isempty, tables));
    if isempty(tables)
        combined = struct();
    elseif isscalar(tables)
        combined = tables{1};
    else
        combined = synchronize(tables{:});
    end
end

function difference = timetable_difference(actual, estimate, name)
    difference = [];
    if isempty(actual) || isempty(estimate)
        return;
    end

    synced = synchronize(actual, estimate, estimate.Time, 'linear');
    values = table2array(synced(:, 1)) - ...
        table2array(synced(:, 2));
    difference = timetable(synced.Time, values, ...
        'VariableNames', {char(name)});
end
