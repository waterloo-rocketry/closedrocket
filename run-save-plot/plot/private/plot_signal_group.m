function plot_signal_group(ax, signals, bases, titleText, yText, tlim)
%PLOT_SIGNAL_GROUP Plot related scalar or vector timetable signals on one axis.

    hold(ax, "on");

    plotted = false;
    for i = 1:numel(signals)
        signal = signals{i};
        if isempty(signal)
            continue;
        end

        [time, values] = signal_data(signal);
        labels = label_for_series(bases(i), size(values, 2));
        style = "-";
        if i > 1
            style = "--";
        end

        for j = 1:size(values, 2)
            plot(ax, time, values(:, j), style, ...
                "DisplayName", labels(j), "LineWidth", 1.0);
            plotted = true;
        end
    end

    title(ax, titleText, "FontWeight", "normal");
    xlabel(ax, "Time (secs)");
    ylabel(ax, yText);
    grid(ax, "on");
    box(ax, "on");
    ax.XMinorTick = "on";
    ax.YMinorTick = "on";
    ax.LineWidth = 0.8;

    if ~isempty(tlim)
        xlim(ax, tlim(1:2));
    end

    if plotted
        legend(ax, "show", "Location", "northoutside", ...
            "Orientation", "horizontal", "NumColumns", 8);
    end
end
