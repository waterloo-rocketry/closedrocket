function data = load_log_data(source, callerName)
%LOAD_LOG_DATA Load a closedrocket log file or normalize a log struct.

    if nargin < 2
        callerName = "plot";
    end

    if isstring(source) || ischar(source)
        loaded = load(source);
        data.file = string(source);
        if isfield(loaded, "log")
            data.log = loaded.log;
            data.signals = loaded.log.signals;
            data.time = getLogTime(loaded.log);
        else
            error(callerName + ":NoLogData", ...
                "File '%s' does not contain a log struct.", source);
        end
    elseif isstruct(source)
        data.file = "";
        data.log = source;
        if isfield(source, "signals")
            data.signals = source.signals;
            data.time = getLogTime(source);
        else
            data.signals = source;
            data.time = [];
        end
    else
        error(callerName + ":UnsupportedInput", ...
            "Expected a saved log file, log struct, or signal struct.");
    end
end

function time = getLogTime(log)
    if isfield(log, "time")
        time = log.time;
    else
        time = [];
    end
end
