function figs = plot_log(logFile, varargin)
%PLOT_LOG Plot a saved simulation log.
%
%   plot_log plots results/sil/result.mat.
%   plot_log(logFile, "PlotError", true) also plots derived nav errors.

    scriptFolder = fileparts(mfilename("fullpath"));

    if nargin < 1 || isempty(logFile)
        folderpath = fullfile(fileparts(scriptFolder), "results", "sil");
        logFile = fullfile(folderpath, "result.mat");
    else
        folderpath = fileparts(char(logFile));
    end

    if isempty(folderpath)
        folderpath = ".";
    end

    parser = inputParser;
    addParameter(parser, "TimeLimits", []);
    addParameter(parser, "PlotError", false);
    addParameter(parser, "SaveFigures", false);
    addParameter(parser, "OutputDir", string(folderpath));
    parse(parser, varargin{:});

    navFigs = plot_navigation(logFile, ...
        "TimeLimits", parser.Results.TimeLimits, ...
        "PlotError", parser.Results.PlotError, ...
        "SaveFigures", parser.Results.SaveFigures, ...
        "OutputDir", parser.Results.OutputDir);

    ctrlFigs = plot_controller(logFile, ...
        "TimeLimits", parser.Results.TimeLimits, ...
        "SaveFigures", parser.Results.SaveFigures, ...
        "OutputDir", parser.Results.OutputDir);

    figs = mergeStructs(navFigs, ctrlFigs);
end

function out = mergeStructs(varargin)
    out = struct();
    for i = 1:nargin
        names = fieldnames(varargin{i});
        for j = 1:numel(names)
            out.(names{j}) = varargin{i}.(names{j});
        end
    end
end
