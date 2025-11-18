function plot_stats_covariance(sdt_array, type, percentiles, varargin)
    % Plot the mean and standard deviation across multiple simulations
    % Input: sdt_array - cell array of sdt structs

    %% Colors (same as plot_state)
    colors(1,:) = [0, 0, 0];      % Black
    colors(2,:) = [0.9, 0.3, 0.1];% Red
    colors(3,:) = [0.2, 0.8, 0.1];% Green
    colors(4,:) = [0.1, 0.2, 0.8];% Blue

    %% Extract valid datasets and reference time
    N = numel(sdt_array);
    valid_idx = [];
    T_ref = [];
    tlim = [0 0];

    for n = 1:N
        if ~isfield(sdt_array{n}, type) || isempty(sdt_array{n}.(type))
            continue
        end
        T_now = seconds(sdt_array{n}.(type).Time);
        if isempty(T_now), continue; end
        valid_idx(end+1) = n;

        if isempty(T_ref) || length(T_now) > length(T_ref)
            T_ref = T_now;
        end
        tlim(2) = max(tlim(2), T_now(end));
    end

    if isempty(valid_idx)
        warning("No valid simulations for type '%s'.", type);
        return;
    end


    %% Override time limits if provided
    if nargin >= 4 && ~isempty(varargin{1})
        tlim = varargin{1};
    end
    
    %% fields and dimensions
    fields = {'P_norm'};
    names = {"det(P)", "cond(P)", "norm(\textit{\textbf{P}})"};
    dims = 3;

  
    %% Initialize storage
    for i = 1:numel(fields)
        all_data.(fields{i}) = nan(length(T_ref), dims(i), numel(valid_idx));
    end
    
    %% Gather data
    for i = 1:numel(valid_idx)
        k = valid_idx(i);
        data = sdt_array{k}.(type);
        for i = 1:numel(fields)
            fld = fields{i};
            D = data.(fld);
            all_data.(fld)(1:length(D),:,k) = D;
        end
    end


    %% Plot selection
    if nargin >= 5
        idx = varargin{2};
    else
        idx = 1:3;
    end

    %% Tiled layout
    if isscalar(idx)
        tiledlayout(1,1,'TileSpacing','Compact');
    elseif numel(idx) <= 2
        tiledlayout(2,1,'TileSpacing','Compact');
    else
        tiledlayout(3,1,'TileSpacing','Compact');
    end
    
    for i = idx
        switch i
            case 1, plots.det      = nexttile;
            case 2, plots.cond     = nexttile;
            case 3, plots.norm     = nexttile;
        end
    end

    %% Plotting
    for i = idx
        fld = fields{1};
        
        color = colors(4,:);

        skip_every = 10;

        data_every = squeeze(all_data.(fld)(:,i,:));  % [time x runs]
        data = data_every(1:skip_every:end,:);
        time = T_ref(1:skip_every:end);
        
        med = median(data,2,"omitmissing");
        Lmid = prctile(data,(100-percentiles(1))/2, 2);
        Umid = prctile(data,100-(100-percentiles(1))/2, 2);
        Lout = prctile(data,(100-percentiles(2))/2, 2);
        Uout = prctile(data,100-(100-percentiles(2))/2, 2);
        
        plotnames = {'det', 'cond', 'norm'};
        ax = plots.(plotnames{i});
        hold(ax, 'on');

        % Outer envelope
        fill(ax, [time; flipud(time)], [Uout; flipud(Lout)], ...
            color, 'FaceAlpha', 0.2, 'EdgeColor', 'none','HandleVisibility','off');
        % Inner envelope
        if percentiles(1) > 0
            plot(ax, time, Umid, ':', 'Color', [color, 0.8],'HandleVisibility','off');
            plot(ax, time, Lmid, ':', 'Color', [color, 0.8],'HandleVisibility','off');
        end

        % Median line
        plot(ax, time, med, 'Color', [color, 0.5], 'LineWidth', 1, 'DisplayName', names{i});
        
        %% Labels & formatting
        title(ax, names{i}, 'FontWeight','Normal');
        % legend(ax, 'show', 'Orientation','vertical', Location='eastoutside', IconColumnWidth=12 )
        grid(ax,'on');
        xlim(ax, tlim);
        xlabel(ax, "Time [s]");
        ax.YAxis.Exponent = 0;
    end
end
