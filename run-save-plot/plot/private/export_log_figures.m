function export_log_figures(figs, outputDir)
%EXPORT_LOG_FIGURES Export each figure in a figure struct as a PDF.

    outputDir = string(outputDir);
    if ~exist(outputDir, "dir")
        mkdir(outputDir);
    end

    names = fieldnames(figs);
    for i = 1:numel(names)
        exportgraphics(figs.(names{i}), ...
            fullfile(outputDir, "log_" + names{i} + ".pdf"), ...
            "ContentType", "vector");
    end
end
