function baseline = save_baseline(outputFile)
%SAVE_BASELINE Configure the plant model and save its baseline variables.
%
%   baseline = save_baseline
%   baseline = save_baseline("results/sil/plant_model_baseline.mat")

    scriptFolder = fileparts(mfilename("fullpath"));
    repoRoot = fileparts(scriptFolder);

    if nargin < 1 || isempty(outputFile)
        outputFile = fullfile(repoRoot, "results", "sil", ...
            "plant_model_baseline.mat");
    end

    outputFile = string(outputFile);
    outputFolder = fileparts(outputFile);
    if strlength(outputFolder) > 0 && ~exist(outputFolder, "dir")
        mkdir(outputFolder);
    end

    oldFolder = pwd;
    cleanup = onCleanup(@() cd(oldFolder));
    cd(repoRoot);

    variablesBefore = string(who);
    configure_plant_model;
    variablesAfter = string(who);
    baselineVariables = setdiff(variablesAfter, ...
        [variablesBefore; "variablesBefore"], "stable");

    if isempty(baselineVariables)
        error("save_baseline:NoVariables", ...
            "The plant configuration script did not create any variables.");
    end

    save(outputFile, baselineVariables{:});
    baseline = load(outputFile);
end
