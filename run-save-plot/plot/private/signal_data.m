function [time, values] = signal_data(signal)
%SIGNAL_DATA Return seconds and 2-D data from a single-variable timetable.

    name = signal.Properties.VariableNames{1};
    time = seconds(signal.Time);
    values = signal.(name);
    values = reshape(values, size(values, 1), []);
end
