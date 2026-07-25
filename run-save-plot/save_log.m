function log = save_log(source, outputFile, sampleRateHz, varargin)
%SAVE_LOG Save every logsout signal at a fixed sample rate.
%
%   log = save_log(out, "result.mat", 100)
%   log = save_log("out.mat", "result.mat", 100)
%
%   The saved file contains one struct, log. log.time stores the shared time
%   vector in seconds, log.signals stores numeric signal arrays, and
%   log.parameterOverrides stores explicit per-run parameter changes.

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
    addParameter(parser, "ParameterOverrides", struct(), ...
        @(value) isstruct(value) && isscalar(value));
    addParameter(parser, "ParameterFile", "");
    addParameter(parser, "Mode", "");
    addParameter(parser, "ModelName", "");
    parse(parser, varargin{:});

    if ~isnumeric(sampleRateHz) || ~isscalar(sampleRateHz) || ...
            ~isfinite(sampleRateHz) || sampleRateHz <= 0
        error("save_log:InvalidSampleRate", ...
            "Sample rate must be a positive numeric scalar in hertz.");
    end

    [outputFile, outputFolder, fileName] = resolveOutputFile(outputFile);

    source = loadSourceIfNeeded(source);
    if isa(source, "Simulink.SimulationOutput") && ~isempty(source.ErrorMessage)
        error("save_log:SimulationFailed", ...
            "Simulation error: %s", source.ErrorMessage);
    end

    logsout = getLogsout(source);
    [signalNames, signalValues] = flattenLogsout(logsout);

    signalTables = cell(numel(signalNames), 1);
    validNames = strings(numel(signalNames), 1);
    validFields = strings(numel(signalNames), 1);
    validCount = 0;
    for i = 1:numel(signalNames)
        signalName = signalNames(i);
        try
            field = signalField(signalName);
            if any(validFields(1:validCount) == field)
                warning("save_log:DuplicateSignal", ...
                    "Skipping duplicate logged signal '%s'.", signalName);
                continue;
            end

            validCount = validCount + 1;
            validNames(validCount) = signalName;
            validFields(validCount) = field;
            signalTables{validCount} = valueToTimetable( ...
                signalValues{i}, signalName);
        catch err
            warning("save_log:SkippedSignal", ...
                "Skipping logged signal '%s': %s", signalName, err.message);
        end
    end

    if validCount == 0
        error("save_log:NoSignals", "No named logged signals could be saved.");
    end

    validNames = validNames(1:validCount);
    validFields = validFields(1:validCount);
    signalTables = signalTables(1:validCount);
    stopTime = simulationStopTime(source, signalTables);
    time = sharedTimeVector(signalTables, stopTime, sampleRateHz);

    signals = struct();
    savedNames = strings(validCount, 1);
    savedCount = 0;
    for i = 1:validCount
        try
            values = timetableDataAtTimes(signalTables{i}, seconds(time));
            signals.(validFields(i)) = values;
            savedCount = savedCount + 1;
            savedNames(savedCount) = validNames(i);
        catch err
            warning("save_log:SkippedSignal", ...
                "Skipping logged signal '%s': %s", validNames(i), err.message);
        end
    end
    savedNames = savedNames(1:savedCount);
    if savedCount == 0
        error("save_log:NoSignals", ...
            "No logged signals could be resampled and saved.");
    end

    log = struct();
    log.format = "closedrocket-log";
    log.formatVersion = 1;
    log.createdOn = datetime("now");
    log.filePath = workspaceRelativePath(outputFolder);
    log.fileName = fileName;
    log.sampleRateHz = sampleRateHz;
    log.stopTime = stopTime;
    log.time = time;
    log.signalNames = savedNames;
    log.signals = signals;
    log.mode = string(parser.Results.Mode);
    log.modelName = string(parser.Results.ModelName);
    log.parameterFile = workspaceRelativeFile(parser.Results.ParameterFile);
    log.parameterOverrides = parser.Results.ParameterOverrides;

    if ~exist(outputFolder, "dir")
        mkdir(outputFolder);
    end

    save(outputFile, "log", "-v7.3");
end

function [outputFile, outputFolder, fileName] = resolveOutputFile(outputFile)
    outputFile = string(outputFile);
    if ~isscalar(outputFile) || ismissing(outputFile) || strlength(outputFile) == 0
        error("save_log:InvalidOutputFile", ...
            "Output file must be a nonempty scalar path.");
    end

    pathParts = split(replace(outputFile, "\", "/"), "/");
    if any(pathParts == "..")
        error("save_log:InvalidOutputFile", ...
            "Output file must not contain parent-folder references.");
    end

    isAbsolute = startsWith(outputFile, filesep) || ...
        startsWith(outputFile, "\\") || ...
        ~isempty(regexp(outputFile, "^[A-Za-z]:[\\/]", "once"));
    if ~isAbsolute
        outputFile = string(fullfile(pwd, outputFile));
    end

    [outputFolder, name, extension] = fileparts(outputFile);
    fileName = name + extension;
end

function relativePath = workspaceRelativePath(outputFolder)
    workspaceRoot = string(fileparts(fileparts(mfilename("fullpath"))));
    workspacePrefix = workspaceRoot + filesep;

    if strcmpi(outputFolder, workspaceRoot)
        relativePath = ".";
    elseif startsWith(outputFolder, workspacePrefix, "IgnoreCase", true)
        relativePath = extractAfter(outputFolder, strlength(workspacePrefix));
    else
        error("save_log:OutputOutsideWorkspace", ...
            "Output file must be inside the workspace '%s'.", workspaceRoot);
    end
end

function relativeFile = workspaceRelativeFile(filePath)
    filePath = string(filePath);
    if ~isscalar(filePath) || ismissing(filePath)
        error("save_log:InvalidParameterFile", ...
            "ParameterFile must be a scalar path.");
    end
    if strlength(filePath) == 0
        relativeFile = "";
        return;
    end

    pathParts = split(replace(filePath, "\", "/"), "/");
    if any(pathParts == "..")
        error("save_log:InvalidParameterFile", ...
            "ParameterFile must not contain parent-folder references.");
    end

    isAbsolute = startsWith(filePath, filesep) || ...
        startsWith(filePath, "\\") || ...
        ~isempty(regexp(filePath, "^[A-Za-z]:[\\/]", "once"));
    if ~isAbsolute
        relativeFile = filePath;
        return;
    end

    [parameterFolder, name, extension] = fileparts(filePath);
    relativeFolder = workspaceRelativePath(parameterFolder);
    if relativeFolder == "."
        relativeFile = name + extension;
    else
        relativeFile = fullfile(relativeFolder, name + extension);
    end
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

function [names, values] = flattenLogsout(logsout)
    if isa(logsout, "Simulink.SimulationData.Dataset")
        [names, values] = flattenDataset(logsout);
    elseif isstruct(logsout)
        names = string(fieldnames(logsout));
        values = cell(numel(names), 1);
        for i = 1:numel(names)
            values{i} = logsout.(char(names(i)));
        end
    else
        error("save_log:UnsupportedLogsout", ...
            "Unsupported logsout type '%s'.", class(logsout));
    end
end

function [names, values] = flattenDataset(dataset)
    nameChunks = cell(dataset.numElements, 1);
    valueChunks = cell(dataset.numElements, 1);
    nameChunks(:) = {strings(0, 1)};
    valueChunks(:) = {cell(0, 1)};
    elementNames = string(getElementNames(dataset));

    for i = 1:dataset.numElements
        element = getElement(dataset, i);

        if isa(element, "Simulink.SimulationData.Dataset")
            [nestedNames, nestedValues] = flattenDataset(element);
            nameChunks{i} = nestedNames;
            valueChunks{i} = nestedValues;
            continue;
        end

        if ~isprop(element, "Values")
            warning("save_log:UnsupportedElement", ...
                "Skipping logged element %d of type '%s'.", i, class(element));
            continue;
        end

        signalName = "";
        if isprop(element, "Name")
            signalName = string(element.Name);
        end
        if strlength(signalName) == 0 && numel(elementNames) >= i
            signalName = elementNames(i);
        end
        if strlength(signalName) == 0
            continue;
        end

        nameChunks{i} = signalName;
        valueChunks{i} = {element.Values};
    end

    names = vertcat(nameChunks{:});
    values = vertcat(valueChunks{:});
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

function stopTime = simulationStopTime(source, signalTables)
    simulationTime = [];
    if isa(source, "Simulink.SimulationOutput")
        try
            simulationTime = source.tout;
        catch
        end
    elseif isstruct(source) && isfield(source, "tout")
        simulationTime = source.tout;
    end

    if isempty(simulationTime)
        stopTime = max(cellfun( ...
            @(tt) seconds(tt.Time(end)), signalTables));
    else
        if istimetable(simulationTime)
            simulationTime = simulationTime.Time;
        elseif isa(simulationTime, "timeseries")
            simulationTime = simulationTime.Time;
        elseif isstruct(simulationTime) && isfield(simulationTime, "Time")
            simulationTime = simulationTime.Time;
        end
        if isduration(simulationTime)
            simulationTime = seconds(simulationTime);
        end
        stopTime = simulationTime(end);
    end

    if ~isnumeric(stopTime) || ~isscalar(stopTime) || ...
            ~isfinite(stopTime) || stopTime <= 0
        error("save_log:InvalidStopTime", ...
            "Stop time must be a positive numeric scalar in seconds.");
    end
    stopTime = double(stopTime);
end

function time = sharedTimeVector(signalTables, stopTime, sampleRateHz)
    startTime = min(cellfun( ...
        @(tt) seconds(tt.Time(1)), signalTables));
    if stopTime < startTime
        error("save_log:InvalidStopTime", ...
            "Stop time must not precede the first logged sample.");
    end

    sampleCount = floor((stopTime - startTime) * sampleRateHz + 1e-9) + 1;
    time = startTime + (0:sampleCount - 1).' / sampleRateHz;
end

function values = timetableDataAtTimes(tt, newTime)
    try
        tt = retime(tt, newTime, "linear");
    catch
        tt = retime(tt, newTime, "nearest");
    end
    values = table2array(tt);
    values = reshape(values, height(tt), []);
end

function field = signalField(signalName)
    field = string(matlab.lang.makeValidName(char(signalName)));
end
