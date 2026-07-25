%% Settings
mode = "sil";
resultPath = mode;             % Subfolder within results
fileName = "result.mat";
sampleRateHz = 100;
plotError = false;
saveFigures = false;

%% Run
[simout, runInfo] = sim_run("Mode", mode);

%% Save
repoRoot = fileparts(mfilename("fullpath"));
resultFolder = fullfile(repoRoot, "results", resultPath);
logFile = fullfile(resultFolder, fileName);

log = save_log(simout, logFile, sampleRateHz, ...
    "InputVariables", runInfo.inputVariables, ...
    "BaselineFile", runInfo.baselineFile, ...
    "Mode", runInfo.mode, ...
    "ModelName", runInfo.modelName);

%% Plot
figs = plot_log(logFile, ...
    "PlotError", plotError, ...
    "SaveFigures", saveFigures, ...
    "OutputDir", resultFolder);
