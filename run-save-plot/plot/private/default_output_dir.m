function outputDir = default_output_dir(source, outputDir)
%DEFAULT_OUTPUT_DIR Resolve figure export folder.

    outputDir = string(outputDir);
    if outputDir == "" && (isstring(source) || ischar(source))
        outputDir = string(fileparts(char(source)));
    end
    if outputDir == ""
        outputDir = ".";
    end
end
