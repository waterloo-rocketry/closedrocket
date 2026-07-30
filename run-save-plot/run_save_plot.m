function log = run_save_plot(mode, varargin)
%RUN_SAVE_PLOT Run a simulation, save its log, and plot the result.
%
%   run_save_plot
%   run_save_plot("sil", 100)
%   run_save_plot("hil", 120, ...
%       "ResultPath", "hil/test-01", ...
%       "FileName", "flight.mat", ...
%       "PlotMode", "live")
%
%   PlotMode may be "live", "log", "both" (default), or "none".
%   ParameterOverrides is a struct of per-run changes to configured values.
%   Call sim_run, save_log, or plot_log directly when outputs are required.

    if nargin < 1 || isempty(mode)
        mode = "sil";
    end

    parser = inputParser;
    addOptional(parser, "StopTime", 180, @isPositiveScalar);
    addParameter(parser, "ResultPath", mode, @isTextScalar);
    addParameter(parser, "FileName", "result.mat", @isTextScalar);
    addParameter(parser, "PlotMode", "both", @isTextScalar);
    addParameter(parser, "ParameterOverrides", struct(), ...
        @(value) isstruct(value) && isscalar(value));
    parse(parser, varargin{:});

    stopTime = parser.Results.StopTime;
    resultPath = string(parser.Results.ResultPath);
    fileName = string(parser.Results.FileName);
    plotMode = lower(string(parser.Results.PlotMode));
    validPlotModes = ["live", "log", "both", "none"];
    if ~ismember(plotMode, validPlotModes)
        error("run_save_plot:InvalidPlotMode", ...
            "PlotMode must be 'live', 'log', 'both', or 'none'.");
    end
    showLivePlots = ismember(plotMode, ["live", "both"]);
    plotSavedLog = ismember(plotMode, ["log", "both"]);

    sampleRateHz = 100;
    plotError = false;
    saveFigures = false;

    %% Run
    [simout, runInfo] = sim_run( ...
        "Mode", mode, ...
        "StopTime", stopTime, ...
        "LivePlots", showLivePlots, ...
        "ParameterOverrides", parser.Results.ParameterOverrides);

    %% Save
    repoRoot = fileparts(fileparts(mfilename("fullpath")));
    resultFolder = fullfile(repoRoot, "results", resultPath);
    logFile = fullfile(resultFolder, fileName);

    log = save_log(simout, logFile, sampleRateHz, ...
        "ParameterOverrides", runInfo.parameterOverrides, ...
        "ParameterFile", runInfo.parameterFile, ...
        "Mode", runInfo.mode, ...
        "ModelName", runInfo.modelName);

    %% Plot
    if plotSavedLog
        plot_log(logFile, ...
            "PlotError", plotError, ...
            "SaveFigures", saveFigures, ...
            "OutputDir", resultFolder);
    end
end

function valid = isPositiveScalar(value)
    valid = isnumeric(value) && isscalar(value) && ...
        isfinite(value) && value > 0;
end

function valid = isTextScalar(value)
    valid = (isstring(value) && isscalar(value)) || ...
        (ischar(value) && isrow(value));
end
