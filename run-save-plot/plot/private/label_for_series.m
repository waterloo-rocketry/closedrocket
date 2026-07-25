function labels = label_for_series(base, width)
%LABEL_FOR_SERIES Build Simulink-scope-like legend labels.

    labels = strings(1, width);
    if width == 1
        labels(1) = base;
    else
        for i = 1:width
            labels(i) = base + ":" + i;
        end
    end
end
