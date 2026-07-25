function signal = first_log_column(signal)
%FIRST_LOG_COLUMN Keep only the first column of a vector-valued signal.

    if isempty(signal)
        return;
    end

    [time, values] = signal_data(signal);
    if size(values, 2) > 1
        signal = struct("Time", time, "Values", values(:, 1));
    end
end
