function labels = label_for_series(base, width)
%LABEL_FOR_SERIES Build component-based legend labels.

    if width == 1
        labels = string(base);
        return;
    end

    switch width
        case 2
            components = ["x", "y"];
        case 3
            components = ["x", "y", "z"];
        case 4
            components = ["w", "x", "y", "z"];
        otherwise
            components = string(1:width);
    end

    labels = string(base) + "_" + components;
end
