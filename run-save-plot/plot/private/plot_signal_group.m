function plot_signal_group(ax, signals, bases, titleText, yText, tlim, colors)
%PLOT_SIGNAL_GROUP Plot related scalar or vector timetable signals on one axis.

    if nargin < 7
        colors = {};
    end

    hold(ax, "on");

    plotted = false;
    for i = 1:numel(signals)
        signal = signals{i};
        if isempty(signal)
            continue;
        end

        [time, values] = signal_data(signal);
        labels = label_for_series(bases(i), size(values, 2));

        for j = 1:size(values, 2)
            line = plot(ax, time, values(:, j), "-", ...
                "DisplayName", labels(j), "LineWidth", 1.5);
            if numel(colors) >= i && size(colors{i}, 1) >= j
                line.Color = colors{i}(j, :);
            end
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
        plotLegend = legend(ax, "show", ...
            "Location", "best", ...
            "NumColumns", 2, ...
            "Interpreter", "none");
        plotLegend.ItemTokenSize = [12, 18];
    end
end
