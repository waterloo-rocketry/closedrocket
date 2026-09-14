function [plots] = plot_stats_state(sdt_array, type, percentiles, varargin)
% Plot statistical envelopes (median, percentile bounds) across multiple datasets
% Follows the same structure and plotting style as the legacy plot_state_old()

    plots = struct();

    %% Colors 
    colors(1,:) = [0, 0, 0];      % Black
    colors(2,:) = [0.9, 0.3, 0.1];% Red
    colors(3,:) = [0.2, 0.8, 0.1];% Green
    colors(4,:) = [0.1, 0.2, 0.8];% Blue

    %% Extract valid datasets and reference time
    N = numel(sdt_array);
    valid_idx = zeros(1, N);
    valid_count = 0;
    T_ref = [];
    tlim = [0 0];

    for n = 1:N
        sdt = sdt_array{n};
        if ~isstruct(sdt) || ~isfield(sdt, type) || ...
                ~istimetable(sdt.(type)) || isempty(sdt.(type))
            continue
        end
        T_now = seconds(sdt.(type).Time);
        if isempty(T_now), continue; end
        valid_count = valid_count + 1;
        valid_idx(valid_count) = n;

        if isempty(T_ref) || length(T_now) > length(T_ref)
            T_ref = T_now;
        end
        tlim(2) = max(tlim(2), T_now(end));
    end
    valid_idx = valid_idx(1:valid_count);

    if isempty(valid_idx)
        warning("No valid simulations for type '%s'.", type);
        return;
    end

    %% Override time limits if provided
    if nargin >= 4 && ~isempty(varargin{1})
        tlim = varargin{1};
    end

    %% Fields & dimensions
    fields = {'q', 'w', 'v', 'alt', 'cl', 'delta'};
    variables = {["$q_w$","$q_x$","$q_y$","$q_z$"], ["$\omega_x$","$\omega_y$","$\omega_z$"], ["$v_x$","$v_y$","$v_z$"], "$r_x$", "$C_{L_\delta}$", "$\delta$"}; 
    names = {'Quaternion', 'Angular velocity [rad/s]', 'Velocity [m/s]', ...
             'Altitude [km]', 'Canard Coefficient', 'Canard Angle [deg]'};
    dims = [4,3,3,1,1,1];

    %% Create storage
    for i = 1:numel(fields)
        all_data.(fields{i}) = nan(length(T_ref), dims(i), numel(valid_idx));
    end

    %% Collect data
    for run_pos = 1:numel(valid_idx)
        k = valid_idx(run_pos);
        data = sdt_array{k}.(type);
        for field_pos = 1:numel(fields)
            fld = fields{field_pos};
            if ~any(strcmp(data.Properties.VariableNames, fld))
                continue;
            end
            D = data.(fld);
            %%% Special unit conversions
            if strcmp(fld,'alt')
                % Convert m → km
               D = D / 1000;
            end
            if strcmp(fld,'delta')
                % rad → deg
                D = rad2deg(D);
            end
            all_data.(fld)(1:length(D),:,run_pos) = D;
        end
    end
    

    %% Plot selection
    if nargin >= 5
        idx = varargin{2};
    else
        idx = 1:6;
    end

    %% Tiled layout
    if isscalar(idx)
        tiledlayout(1,1,'TileSpacing','Compact');
    elseif numel(idx) <= 2
        tiledlayout(1,2,'TileSpacing','Compact');
    elseif numel(idx) <= 4
        tiledlayout(2,2,'TileSpacing','Compact');
    else
        tiledlayout(3,2,'TileSpacing','Compact');
    end

    %% Tile mapping
    for i = idx
        switch i
            case 1, plots.q     = nexttile;
            case 2, plots.w     = nexttile;
            case 3, plots.v     = nexttile;
            case 4, plots.alt   = nexttile;
            case 5, plots.cl    = nexttile;
            case 6, plots.delta = nexttile;
        end


        %% Plot for each selected field
        fld = fields{i};
        ax = plots.(fld);
        hold(ax, 'on');
        dim = dims(i);

        switch dim
            case 1, color = colors(4,:);
            case 2, color = colors(3:4,:);
            case 3, color = colors(2:4,:);
            case 4, color = colors(1:4,:);
        end

        for  d = 1:dim % 1:dims_plot

            skip_every = 10;

            data_every = squeeze(all_data.(fld)(:,d,:)); % [time × runs]
            data = data_every(1:skip_every:end,:);
            time = T_ref(1:skip_every:end);

            med = median(data,2,"omitmissing");
            Lmid = prctile(data,(100-percentiles(1))/2, 2);
            Umid = prctile(data,100-(100-percentiles(1))/2, 2);
            Lout = prctile(data,(100-percentiles(2))/2, 2);
            Uout = prctile(data,100-(100-percentiles(2))/2, 2);

            % Outer envelope
            fill(ax, [time; flipud(time)], [Uout; flipud(Lout)], ...
                color(d,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off');
            % Inner envelope
            if percentiles(1) > 0
                plot(ax, time, Umid, ':', 'Color', [color(d,:), 0.8],'HandleVisibility','off');
                plot(ax, time, Lmid, ':', 'Color', [color(d,:), 0.8],'HandleVisibility','off');
            end

            % Median line
            plot(ax, time, med, 'Color', [color(d,:), 0.5], 'LineWidth', 1, 'DisplayName', variables{i}(d));
        end
        hold(ax, 'off');

        %% Labels & formatting
        title(ax, names{i}, 'FontWeight','Normal');
        legend(ax, 'show', 'Orientation','vertical', Location='eastoutside', IconColumnWidth=12 )
        grid(ax,'on');
        xlim(ax, tlim);
        % xticks(ax, linspace(tlim(1), tlim(2), 6))
        xlabel(ax, "Time [s]");
        ax.YAxis.Exponent = 0;
        %%% only plot some values
        % if dims_plot == 3
        %     ylim(ax, [-2, 2])
        % end
        % dims_plot = dims_plot +2
    end
end
