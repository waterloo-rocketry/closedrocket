function signal = get_log_signal(data, varargin)
%GET_LOG_SIGNAL Return the first existing signal from log.signals.

    signal = [];
    for i = 1:numel(varargin)
        field = matlab.lang.makeValidName(char(varargin{i}));
        if isfield(data.signals, field)
            storedSignal = data.signals.(field);
            if isnumeric(storedSignal) && ~isempty(data.time)
                signal = struct("Time", data.time, "Values", storedSignal);
            else
                signal = storedSignal;
            end
        end
        if ~isempty(signal)
            return;
        end
    end
end
