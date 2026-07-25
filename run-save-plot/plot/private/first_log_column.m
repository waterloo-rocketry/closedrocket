function signal = first_log_column(signal)
%FIRST_LOG_COLUMN Keep only the first column of a vector-valued signal.

    if isempty(signal)
        return;
    end

    name = signal.Properties.VariableNames{1};
    values = signal.(name);
    if size(values, 2) > 1
        signal = timetable(signal.Time, values(:, 1), 'VariableNames', {name});
    end
end
