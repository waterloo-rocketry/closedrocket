function [simout, runInfo] = sim_run(varargin)
%SIM_RUN Configure and run a SIL or HIL simulation.
%
%   [simout, runInfo] = sim_run
%   [simout, runInfo] = sim_run("StopTime", 120)
%   [simout, runInfo] = sim_run("Mode", "hil", "StopTime", 120, ...
%       "LivePlots", false)
%   [simout, runInfo] = sim_run("ParameterOverrides", ...
%       struct("wind_const_strength", 15))
%
%   runInfo contains the model and input metadata used by save_log.

    modeParser = inputParser;
    modeParser.KeepUnmatched = true;
    addParameter(modeParser, "Mode", "sil");
    parse(modeParser, varargin{:});

    mode = lower(string(modeParser.Results.Mode));
    repoRoot = fileparts(fileparts(mfilename("fullpath")));
    defaults = modeDefaults(mode, repoRoot);

    parser = inputParser;
    addParameter(parser, "Mode", mode);
    addParameter(parser, "StopTime", 180);
    addParameter(parser, "ShowProgress", true);
    addParameter(parser, "LivePlots", true, ...
        @(value) islogical(value) && isscalar(value));
    addParameter(parser, "ParameterOverrides", struct(), ...
        @(value) isstruct(value) && isscalar(value));
    parse(parser, varargin{:});

    stopTime = parser.Results.StopTime;
    if ~isnumeric(stopTime) || ~isscalar(stopTime) || ...
            ~isfinite(stopTime) || stopTime <= 0
        error("sim_run:StopTimeRequired", ...
            "StopTime must be a positive numeric scalar in seconds.");
    end

    oldFolder = pwd;
    cleanup = onCleanup(@() cd(oldFolder));
    cd(repoRoot);

    parameterFile = defaults.parameterFile;
    parameterFolder = fileparts(parameterFile);
    if ~exist(parameterFolder, "dir")
        mkdir(parameterFolder);
    end

    save_sim_parameters(parameterFile);

    simin = Simulink.SimulationInput(defaults.modelName);
    simin = setModelParameter(simin, "StopTime", string(stopTime));
    simin = simin.loadVariablesFromMATFile(parameterFile);
    simin = applyParameterOverrides( ...
        simin, parser.Results.ParameterOverrides);

    if ~parser.Results.LivePlots
        scopeCleanup = suppressLiveScopes(defaults.modelName); %#ok<NASGU>
    end

    callbackVariables = ["elems", "Sensor1DBus", "Sensor3DBus"];
    callbackState = captureBaseVariables(callbackVariables);
    callbackCleanup = onCleanup(@() restoreBaseVariables( ...
        callbackState));

    progress = "off";
    if parser.Results.ShowProgress
        progress = "on";
    end

    simout = sim(simin, "ShowProgress", progress);

    runInfo = struct();
    runInfo.mode = mode;
    runInfo.modelName = defaults.modelName;
    runInfo.stopTime = stopTime;
    runInfo.livePlots = logical(parser.Results.LivePlots);
    runInfo.parameterFile = string(parameterFile);
    runInfo.defaultLogFile = string(defaults.logFile);
    runInfo.parameterOverrides = parser.Results.ParameterOverrides;
end

function defaults = modeDefaults(mode, repoRoot)
    switch mode
        case "sil"
            defaults.modelName = "plant-model/CC_Flight_Simulation";
            defaults.logFile = fullfile(repoRoot, "results", "sil", ...
                "result.mat");
            defaults.parameterFile = fullfile( ...
                "results", "sil", "sim_parameters.mat");
        case "hil"
            defaults.modelName = "hil/hil";
            defaults.logFile = fullfile(repoRoot, "results", "hil", ...
                "result.mat");
            defaults.parameterFile = fullfile( ...
                "results", "hil", "sim_parameters.mat");
        otherwise
            error("sim_run:UnknownMode", ...
                "Mode must be 'sil' or 'hil', not '%s'.", mode);
    end
end

function simin = applyParameterOverrides(simin, overrides)
    if isempty(overrides)
        return;
    end

    names = fieldnames(overrides);
    for i = 1:numel(names)
        simin = simin.setVariable(names{i}, overrides.(names{i}));
    end
end

function cleanup = suppressLiveScopes(modelFile)
    modelHandle = load_system(modelFile);
    modelName = getfullname(modelHandle);
    scopeBlocks = find_system(modelName, ...
        "LookUnderMasks", "all", ...
        "FollowLinks", "on", ...
        "BlockType", "Scope");
    openStates = get_param(scopeBlocks, "OpenAtSimulationStart");
    liveScopes = scopeBlocks(strcmp(openStates, "on"));
    liveStates = openStates(strcmp(openStates, "on"));
    diagrams = find_system("Type", "block_diagram");
    dirtyStates = get_param(diagrams, "Dirty");

    for i = 1:numel(liveScopes)
        set_param(liveScopes{i}, "OpenAtSimulationStart", "off");
        set_param(liveScopes{i}, "Open", "off");
    end

    dirtyAfter = get_param(diagrams, "Dirty");
    changed = ~strcmp(dirtyStates, dirtyAfter);
    changedDiagrams = diagrams(changed);
    changedDirtyStates = dirtyStates(changed);
    restoreDirtyStates(changedDiagrams, changedDirtyStates);

    cleanup = onCleanup(@() restoreLiveScopes( ...
        liveScopes, liveStates, changedDiagrams, changedDirtyStates));
end

function restoreLiveScopes(scopeBlocks, openStates, diagrams, dirtyStates)
    for i = 1:numel(scopeBlocks)
        if getSimulinkBlockHandle(scopeBlocks{i}) ~= -1
            set_param(scopeBlocks{i}, ...
                "OpenAtSimulationStart", openStates{i});
        end
    end
    restoreDirtyStates(diagrams, dirtyStates);
end

function restoreDirtyStates(diagrams, dirtyStates)
    for i = 1:numel(diagrams)
        if bdIsLoaded(diagrams{i})
            set_param(diagrams{i}, "Dirty", dirtyStates{i});
        end
    end
end

function state = captureBaseVariables(variableNames)
    state = repmat(struct( ...
        "name", "", ...
        "existed", false, ...
        "value", []), numel(variableNames), 1);

    for i = 1:numel(variableNames)
        state(i).name = variableNames(i);
        state(i).existed = evalin("base", ...
            "exist('" + variableNames(i) + "', 'var') == 1");
        if state(i).existed
            state(i).value = evalin("base", variableNames(i));
        end
    end
end

function restoreBaseVariables(state)
    for i = 1:numel(state)
        if state(i).existed
            assignin("base", state(i).name, state(i).value);
        else
            evalin("base", "clear " + state(i).name);
        end
    end
end
