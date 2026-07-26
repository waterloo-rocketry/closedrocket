function err = signal_difference(actual, estimate, variableName)
%SIGNAL_DIFFERENCE Compute actual - estimate on the actual signal timebase.

    err = [];
    if isempty(actual) || isempty(estimate)
        return;
    end

    [actualTime, actualValues] = signal_data(actual);
    [estimateTime, estimateValues] = signal_data(estimate);
    if ~isequal(actualTime, estimateTime)
        estimateValues = interp1(estimateTime, estimateValues, ...
            actualTime, "linear", NaN);
    end
    width = min(size(actualValues, 2), size(estimateValues, 2));

    if width < 1
        return;
    end

    err = struct( ...
        "Time", actualTime, ...
        "Values", actualValues(:, 1:width) - estimateValues(:, 1:width), ...
        "Name", string(variableName));
end
