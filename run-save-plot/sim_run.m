function [simout, runInfo] = sim_run(varargin)
%SIM_RUN Configure and run a SIL or HIL simulation.
%
%   simout = sim_run runs plant-model/CC_Flight_Simulation.
%
%   [simout, runInfo] = sim_run("StopTime", 120)
%   [simout, runInfo] = sim_run("Mode", "hil")
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
    addParameter(parser, "StopTime", defaults.stopTime);
    addParameter(parser, "ShowProgress", true);
    addParameter(parser, "Variables", struct());
    parse(parser, varargin{:});

    oldFolder = pwd;
    cleanup = onCleanup(@() cd(oldFolder));
    cd(repoRoot);

    baselineFile = defaults.baselineFile;
    baselineFolder = fileparts(baselineFile);
    if ~exist(baselineFolder, "dir")
        mkdir(baselineFolder);
    end

    baseline = save_baseline(baselineFile);

    simin = Simulink.SimulationInput(defaults.modelName);
    simin = setModelParameter(simin, "StopTime", string(parser.Results.StopTime));
    simin = simin.loadVariablesFromMATFile(baselineFile);
    simin = applyVariables(simin, parser.Results.Variables);

    progress = "off";
    if parser.Results.ShowProgress
        progress = "on";
    end

    simout = sim(simin, "ShowProgress", progress);

    runInfo = struct();
    runInfo.mode = mode;
    runInfo.modelName = defaults.modelName;
    runInfo.baselineFile = string(baselineFile);
    runInfo.defaultLogFile = string(defaults.logFile);
    runInfo.inputVariables = inputVariableDiff(simin, baseline);
end

function defaults = modeDefaults(mode, repoRoot)
    switch mode
        case "sil"
            defaults.modelName = "plant-model/CC_Flight_Simulation";
            defaults.stopTime = 100;
            defaults.logFile = fullfile(repoRoot, "results", "sil", ...
                "result.mat");
            defaults.baselineFile = fullfile(repoRoot, "results", "sil", ...
                "plant_model_baseline.mat");
        case "hil"
            defaults.modelName = "hil/hil";
            defaults.stopTime = 120;
            defaults.logFile = fullfile(repoRoot, "results", "hil", ...
                "result.mat");
            defaults.baselineFile = fullfile(repoRoot, "results", "hil", ...
                "plant_model_baseline.mat");
        otherwise
            error("sim_run:UnknownMode", ...
                "Mode must be 'sil' or 'hil', not '%s'.", mode);
    end
end

function simin = applyVariables(simin, variables)
    if isempty(variables)
        return;
    end

    names = fieldnames(variables);
    for i = 1:numel(names)
        simin = simin.setVariable(names{i}, variables.(names{i}));
    end
end

function inputVariables = inputVariableDiff(simin, baseline)
    inputVariables = struct();
    for i = 1:length(simin.Variables)
        variable = simin.Variables(i);
        inputVariables.(variable.Name) = variable.Value;
    end

    shared = intersect(fieldnames(inputVariables), fieldnames(baseline));
    for i = 1:numel(shared)
        field = shared{i};
        if isequaln(inputVariables.(field), baseline.(field))
            inputVariables = rmfield(inputVariables, field);
        end
    end
end
