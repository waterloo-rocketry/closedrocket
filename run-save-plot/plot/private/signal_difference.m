function err = signal_difference(actual, estimate, variableName)
%SIGNAL_DIFFERENCE Compute actual - estimate on the actual signal timebase.

    err = [];
    if isempty(actual) || isempty(estimate)
        return;
    end

    [~, actualValues] = signal_data(actual);
    estimate = retime(estimate, actual.Time, "linear");
    [~, estimateValues] = signal_data(estimate);
    width = min(size(actualValues, 2), size(estimateValues, 2));

    if width < 1
        return;
    end

    err = timetable(actual.Time, actualValues(:, 1:width) - estimateValues(:, 1:width), ...
        'VariableNames', cellstr(variableName));
end
