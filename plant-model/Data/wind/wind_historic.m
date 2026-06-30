function [wind_heights, wind_vectors] = wind_historic(wind_layer_threshold)
    %% windspeed wrt height from csv
    wind_data = readtable("plant-model/data/wind/sim_parameters_historical_aug2023+aug2024_no-clouds_2025-06-03.csv", ...
        ReadVariableNames=true, VariableNamingRule="preserve");
    
    %% pick a random timestamp and process
    for iter = 1:100
        row_idx = randi(height(wind_data));
        %row_idx = find(wind_data{:,1}==datetime(2024, 08, 10, 14, 0, 0)); % select specific date
    
        row = wind_data(row_idx, :);
        wind_time = row.date;
        
        cols = wind_data.Properties.VariableNames;
        
        % Speed columns (numeric names, e.g., '5', '10', etc.)
        speed_cols = cols(~cellfun(@isempty, regexp(cols, "\d+$")));
        
        % Direction columns (e.g., 'direction [5m]', 'direction [10m]')
        dir_cols = cols(~cellfun(@isempty, regexp(cols, "direction")));
        
        % Convert column names to numeric heights [m]
        heights = str2double(speed_cols);
        
        % Extract wind speeds and directions for selected timestamp
        speeds = row{1, speed_cols}; % [m/s]
        dirs   = row{1, dir_cols};   % [deg]
    
        if all(speeds <= wind_layer_threshold)
            fprintf('Selected historic wind date: %s\n', wind_time);
            break
        elseif iter == 100
            fprintf('Selected historic wind date: %s\n Speed threshold exceeded after many tries.\n', wind_time);
        else
            % fprintf('try another day\n');
            continue
        end
        
    end
    
    %% Convert to 3D wind vectors
    % im assuming rocket axes (x = up, y = north, z = east)
    % and that the direction is deg from north
    wind_vectors = zeros(length(heights), 3);
    
    for i = 1:length(heights)
        deg = dirs(i);
        rad = deg2rad(deg);
        % stolen from the way wind_const_strength and wind_const_direction were used before
        wind_vectors(i, :) = speeds(i) * [0, cos(rad), sin(rad)]; % no vertical component
    end
    
    % N x 1
    wind_heights = heights';
    % N x 3
    wind_vectors = wind_vectors;
    
    
    % fprintf('heights:\n');
    % disp(wind_heights);
    % fprintf("windvectors:\n");
    % disp(wind_vectors)
end