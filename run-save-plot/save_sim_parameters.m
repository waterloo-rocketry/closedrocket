function parameters = save_sim_parameters(outputFile)
%SAVE_SIM_PARAMETERS Configure the plant model and save its parameters.
%
%   parameters = save_sim_parameters
%   parameters = save_sim_parameters("results/sil/sim_parameters.mat")

    scriptFolder = fileparts(mfilename("fullpath"));
    repoRoot = fileparts(scriptFolder);

    if nargin < 1 || isempty(outputFile)
        outputFile = fullfile(repoRoot, "results", "sil", ...
            "sim_parameters.mat");
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
    parameterVariables = setdiff(variablesAfter, ...
        [variablesBefore; "variablesBefore"], "stable");

    if isempty(parameterVariables)
        error("save_sim_parameters:NoVariables", ...
            "The plant configuration script did not create any variables.");
    end

    save(outputFile, parameterVariables{:});
    if nargout > 0
        parameters = load(outputFile);
    end
end
