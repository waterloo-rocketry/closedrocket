function signal = get_log_signal(data, varargin)
%GET_LOG_SIGNAL Return the first existing timetable from log.signals.

    signal = [];
    for i = 1:numel(varargin)
        field = matlab.lang.makeValidName(char(varargin{i}));
        if isfield(data.signals, field)
            signal = data.signals.(field);
        end
        if ~isempty(signal)
            return;
        end
    end
end
