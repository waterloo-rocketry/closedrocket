function log = save_log(source, outputFile, sampleRateHz, varargin)
%SAVE_LOG Save every logsout signal at a fixed sample rate.
%
%   log = save_log(out, "result.mat", 100)
%   log = save_log("out.mat", "result.mat", 100)
%
%   The saved file contains one struct, log. Each logged signal is stored as
%   a timetable in log.signals using the original logsout signal name.

    if nargin < 1 || isempty(source)
        source = evalin("base", "out");
    end
    if nargin < 2 || isempty(outputFile)
        outputFile = "result.mat";
    end
    if nargin < 3 || isempty(sampleRateHz)
        sampleRateHz = 100;
    end

    parser = inputParser;
    addParameter(parser, "InputVariables", struct());
    addParameter(parser, "BaselineFile", "");
    addParameter(parser, "Mode", "");
    addParameter(parser, "ModelName", "");
    parse(parser, varargin{:});

    source = loadSourceIfNeeded(source);
    if isa(source, "Simulink.SimulationOutput") && ~isempty(source.ErrorMessage)
        error("save_log:SimulationFailed", ...
            "Simulation error: %s", source.ErrorMessage);
    end

    logsout = getLogsout(source);
    signalNames = getLogElementNames(logsout);

    signals = struct();
    for i = 1:numel(signalNames)
        signalName = signalNames(i);
        element = getLogElement(logsout, signalName);
        try
            tt = valueToTimetable(element.Values, signalName);
            signals.(signalField(signalName)) = retimeAtRate(tt, sampleRateHz);
        catch err
            warning("save_log:SkippedSignal", ...
                "Skipping logged signal '%s': %s", signalName, err.message);
        end
    end

    log = struct();
    log.format = "closedrocket-log";
    log.formatVersion = 1;
    log.createdOn = datetime("now");
    log.sampleRateHz = sampleRateHz;
    log.signalNames = signalNames;
    log.signals = signals;
    log.mode = string(parser.Results.Mode);
    log.modelName = string(parser.Results.ModelName);
    log.baselineFile = string(parser.Results.BaselineFile);
    log.inputVariables = parser.Results.InputVariables;

    outputFolder = fileparts(char(outputFile));
    if ~isempty(outputFolder) && ~exist(outputFolder, "dir")
        mkdir(outputFolder);
    end

    save(outputFile, "log", "-v7.3");
end

function source = loadSourceIfNeeded(source)
    if isstring(source) || ischar(source)
        loaded = load(source);
        if isfield(loaded, "out")
            source = loaded.out;
        elseif isfield(loaded, "simout")
            source = loaded.simout;
        elseif isfield(loaded, "logsout")
            source = loaded.logsout;
        elseif isfield(loaded, "log")
            source = loaded.log;
        else
            names = fieldnames(loaded);
            source = loaded.(names{1});
        end
    end
end

function logsout = getLogsout(source)
    if isa(source, "Simulink.SimulationData.Dataset")
        logsout = source;
    elseif isa(source, "Simulink.SimulationOutput")
        logsout = source.logsout;
    elseif isstruct(source) && isfield(source, "logsout")
        logsout = source.logsout;
    elseif isstruct(source) && isfield(source, "signals")
        logsout = source.signals;
    else
        error("save_log:NoLogsout", ...
            "Could not find a logsout dataset in the provided source.");
    end
end

function names = getLogElementNames(logsout)
    if isa(logsout, "Simulink.SimulationData.Dataset")
        try
            names = string(getElementNames(logsout));
        catch
            names = strings(1, logsout.numElements);
            for i = 1:logsout.numElements
                element = getElement(logsout, i);
                names(i) = string(element.Name);
            end
        end
    elseif isstruct(logsout)
        names = string(fieldnames(logsout));
    else
        error("save_log:UnsupportedLogsout", ...
            "Unsupported logsout type '%s'.", class(logsout));
    end
end

function element = getLogElement(logsout, signalName)
    if isa(logsout, "Simulink.SimulationData.Dataset")
        element = getElement(logsout, char(signalName));
    elseif isstruct(logsout)
        element.Values = logsout.(signalField(signalName));
    end
end

function tt = valueToTimetable(values, signalName)
    if istimetable(values)
        tt = values;
    elseif isa(values, "timeseries")
        tt = timeseries2timetable(values);
    elseif isstruct(values) && isfield(values, "Time") && isfield(values, "Data")
        time = values.Time;
        if isnumeric(time)
            time = seconds(time);
        end
        tt = timetable(time(:), values.Data, 'VariableNames', cellstr(signalField(signalName)));
    else
        error("save_log:UnsupportedSignalValue", ...
            "Unsupported logged value type '%s'.", class(values));
    end

    tt = compactToSingleVariable(tt, signalField(signalName));
end

function tt = compactToSingleVariable(tt, variableName)
    if isscalar(tt.Properties.VariableNames)
        tt.Properties.VariableNames = cellstr(variableName);
    else
        data = table2array(tt);
        tt = timetable(tt.Time, data, 'VariableNames', cellstr(variableName));
    end
end

function tt = retimeAtRate(tt, sampleRateHz)
    if isempty(tt) || height(tt) == 0 || sampleRateHz <= 0
        return;
    end

    startTime = tt.Time(1);
    stopTime = tt.Time(end);
    newTime = (startTime:seconds(1 / sampleRateHz):stopTime).';
    if isempty(newTime)
        newTime = startTime;
    end

    try
        tt = retime(tt, newTime, "linear");
    catch
        tt = retime(tt, newTime, "nearest");
    end
end

function field = signalField(signalName)
    field = string(matlab.lang.makeValidName(char(signalName)));
end
