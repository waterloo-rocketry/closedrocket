function [time, values] = signal_data(signal)
%SIGNAL_DATA Return seconds and 2-D values from a log signal.

    if istimetable(signal)
        name = signal.Properties.VariableNames{1};
        time = signal.Time;
        values = signal.(name);
    elseif isstruct(signal) && isfield(signal, "Time") && ...
            isfield(signal, "Values")
        time = signal.Time;
        values = signal.Values;
    else
        error("plot_log:UnsupportedSignal", ...
            "Expected a timetable or a signal with Time and Values.");
    end

    if isduration(time)
        time = seconds(time);
    end
    time = time(:);
    values = reshape(values, size(values, 1), []);
end
