function plots = plot_stats_control(sdt_array, type, commontitle, percentiles)
%PLOT_STATS_CONTROL Plot control envelopes across Monte Carlo simulations.

    fields = {'ref', 'err', 'cmd', 'roll', 'rate', 'delta'};
    names = {'Reference [rad]', 'Roll error [rad]', 'Command [rad]', ...
        'Roll angle [rad]', 'Rates [rad/s]', 'Actuation [rad]'};
    dims = [1, 1, 1, 1, 1, 1];

    [valid_idx, T_ref] = valid_runs(sdt_array, type);
    if isempty(valid_idx)
        warning("No valid simulations found for type '%s'. Skipping plot.", type);
        plots = struct();
        return;
    end

    for f = 1:numel(fields)
        all_data.(fields{f}) = nan(length(T_ref), dims(f), numel(valid_idx));
    end

    for run_pos = 1:numel(valid_idx)
        k = valid_idx(run_pos);
        data = sdt_array{k}.(type);
        for f = 1:numel(fields)
            field = fields{f};
            if ~any(strcmp(data.Properties.VariableNames, field))
                continue;
            end
            values = data.(field);
            all_data.(field)(1:length(values), :, run_pos) = values;
        end
    end

    var_colors = [ ...
        0.00, 0.45, 0.74; ...
        0.00, 0.60, 0.20; ...
        1.00, 0.50, 0.00; ...
        0.60, 0.20, 0.50];

    tlo = tiledlayout(2, 3, "TileSpacing", "Compact", "Padding", "Compact");
    plots = struct();

    for f = 1:numel(fields)
        ax = nexttile;
        plots.(fields{f}) = ax;
        hold(ax, "on");

        for d = 1:dims(f)
            color = var_colors(d, :);
            data = squeeze(all_data.(fields{f})(:, d, :));
            med = median(data, 2, "omitmissing");
            lower_mid = prctile(data, (100 - percentiles(1)) / 2, 2);
            upper_mid = prctile(data, 100 - (100 - percentiles(1)) / 2, 2);
            lower = prctile(data, (100 - percentiles(2)) / 2, 2);
            upper = prctile(data, 100 - (100 - percentiles(2)) / 2, 2);

            fill(ax, [T_ref; flipud(T_ref)], [upper; flipud(lower)], ...
                color, "FaceAlpha", 0.15, "EdgeColor", "none", "HandleVisibility", "off");
            plot(ax, T_ref, lower_mid, ":", "Color", color, "LineWidth", 1, "HandleVisibility", "off");
            plot(ax, T_ref, upper_mid, ":", "Color", color, "LineWidth", 1, "HandleVisibility", "off");
            plot(ax, T_ref, med, "Color", color, "LineWidth", 1.5);
        end

        title(ax, names{f}, "FontWeight", "Normal");
        grid(ax, "on");
        xlabel(ax, "Time [s]");
    end

    title(tlo, commontitle);

    ax = plots.(fields{end});
    h_median = plot(ax, nan, nan, "-", "Color", [0 0 0], "LineWidth", 1.5);
    h_perc1 = plot(ax, nan, nan, ":", "Color", [0 0 0], "LineWidth", 1);
    h_fill = fill(ax, nan, nan, [0 0 0], "FaceAlpha", 0.15, "EdgeColor", "none");
    labels = {'Med.', sprintf('%d%%', percentiles(1)), sprintf('%d%%', percentiles(2))};
    lgd = legend(ax, [h_median, h_perc1, h_fill], labels, ...
        "FontSize", 8, "Orientation", "horizontal", "NumColumns", 3);
    set(lgd, "Units", "normalized");
    lgd.Position(1:2) = [0.56, 0.97];
end

function [valid_idx, T_ref] = valid_runs(sdt_array, type)
    valid_idx = [];
    T_ref = [];

    for n = 1:numel(sdt_array)
        sdt = sdt_array{n};
        if ~isstruct(sdt) || ~isfield(sdt, type) || ...
                ~istimetable(sdt.(type)) || isempty(sdt.(type))
            continue;
        end

        T_now = seconds(sdt.(type).Time);
        if isempty(T_now)
            continue;
        end

        valid_idx(end + 1) = n; %#ok<AGROW>
        if isempty(T_ref) || length(T_now) > length(T_ref)
            T_ref = T_now;
        end
    end
end
