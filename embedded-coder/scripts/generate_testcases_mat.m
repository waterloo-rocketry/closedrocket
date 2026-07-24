function generate_testcases_mat(output_file)
%GENERATE_TESTCASES_MAT Convert controller/navigation logs to testcases.mat.
%
% The saved file contains controller_testcases and navigator_testcases.
% Navigator statuses are varied deterministically to exercise sensor dropouts.

script_dir = fileparts(mfilename('fullpath'));
if nargin < 1
    output_file = fullfile(script_dir, 'testcases.mat');
end

controller_raw = parse_log(fullfile(script_dir, 'controller_log.txt'), 'timestamp');
navigator_raw = parse_log(fullfile(script_dir, 'navigator_log.txt'), 'dt');

controller_testcases = repmat(struct(), size(controller_raw));
for i = 1:numel(controller_raw)
    raw = controller_raw(i);
    controller_testcases(i).time = raw.timestamp;
    controller_testcases(i).dt = raw.dt;
    controller_testcases(i).roll_state = raw.roll_state;
    controller_testcases(i).pdyn = raw.pdyn;
    controller_testcases(i).delta = raw.delta;
    controller_testcases(i).ctrl_mem.coeffs = raw.coeffs;
    controller_testcases(i).ctrl_mem.w = raw.w;
    controller_testcases(i).ctrl_mem.P = raw.P;
    controller_testcases(i).ctrl_mem.delta_lp = raw.delta_lp;
    controller_testcases(i).ctrl_mem.w_dot_lp = raw.w_dot_lp;
end

sensor_names = {'board_accel', 'board_gyro', 'mti_accel', ...
    'mti_gyro', 'ad_accel', 'ad_gyro', 'board_baro', ...
    'board_mag', 'mti_baro', 'mti_mag'};
environment_sensor_names = {'board_baro', 'board_mag', ...
    'mti_baro', 'mti_mag'};
status_dropout_probability = 0.02;
status_random_stream = RandStream('mt19937ar', 'Seed', 20260716);
navigator_testcases = struct([]);

for i = 1:numel(navigator_raw)
    raw = navigator_raw(i);
    testcase.dt = raw.dt;
    testcase.flight_phase = logical(raw.flight_phase);
    testcase.x = raw.x;
    testcase.P = raw.P;

    testcase.bias.board_gyro = raw.bias_board_gyro;
    testcase.bias.mti_gyro = raw.bias_mti_gyro;
    testcase.bias.ad_gyro = raw.bias_ad_gyro;
    testcase.bias.board_mag_earth = raw.bias_board_mag;
    testcase.bias.mti_mag_earth = raw.bias_mti_mag;
    testcase.bias.board_baro = raw.bias_board_baro;
    testcase.bias.mti_baro = raw.bias_mti_baro;

    for j = 1:numel(sensor_names)
        sensor = sensor_names{j};
        testcase.sens_filt.(sensor) = raw.([sensor '_filt']);
        testcase.sens_input.(sensor).meas = raw.([sensor '_meas']);
        testcase.sens_input.(sensor).status = logical(raw.([sensor '_status']));
    end

    % Enable all barometers and magnetometers in every second testcase.
    if mod(i, 2) == 0
        for j = 1:numel(environment_sensor_names)
            sensor = environment_sensor_names{j};
            testcase.sens_input.(sensor).status = true;
        end
    end

    % Give every enabled sensor an independent 2%% chance of dropping out.
    for j = 1:numel(sensor_names)
        sensor = sensor_names{j};
        if testcase.sens_input.(sensor).status && ...
                rand(status_random_stream) < status_dropout_probability
            testcase.sens_input.(sensor).status = false;
        end
    end

    navigator_testcases = append_record(navigator_testcases, testcase);
    clear testcase
end

save(output_file, 'controller_testcases', 'navigator_testcases');
fprintf('Saved %d controller and %d navigator testcases to %s\n', ...
    numel(controller_testcases), numel(navigator_testcases), output_file);
end


function records = parse_log(file_path, first_field)
lines = regexp(fileread(file_path), '\r?\n', 'split');
field_pattern = '^([A-Za-z][A-Za-z0-9_]*)\s*=\s*$';
records = struct([]);
record = struct();
field_name = '';
field_lines = {};

for i = 1:numel(lines)
    match = regexp(strtrim(lines{i}), field_pattern, 'tokens', 'once');
    if ~isempty(match)
        [record, field_name, field_lines] = finish_field( ...
            record, field_name, field_lines, file_path);
        next_field = match{1};
        if strcmp(next_field, first_field) && ~isempty(fieldnames(record))
            records = append_record(records, record);
            record = struct();
        end
        field_name = next_field;
    elseif ~isempty(field_name)
        field_lines{end + 1} = strtrim(lines{i}); %#ok<AGROW>
    end
end

[record, ~, ~] = finish_field(record, field_name, field_lines, file_path);
if ~isempty(fieldnames(record))
    records = append_record(records, record);
end
assert(~isempty(records), 'No records beginning with "%s" found in %s.', ...
    first_field, file_path);
end


function records = append_record(records, record)
if isempty(records)
    records = record;
else
    expected_fields = fieldnames(records);
    actual_fields = fieldnames(record);
    assert(isequal(sort(expected_fields), sort(actual_fields)), ...
        'Log records do not all contain the same fields.');
    records(end + 1) = orderfields(record, records); %#ok<AGROW>
end
end


function [record, field_name, field_lines] = finish_field( ...
    record, field_name, field_lines, file_path)
if ~isempty(field_name)
    try
        record.(field_name) = parse_numeric_value(field_lines);
    catch exception
        error('Invalid field "%s" in %s: %s', ...
            field_name, file_path, exception.message);
    end
end
field_name = '';
field_lines = {};
end


function value = parse_numeric_value(lines)
% MATLAB displays wide matrices in successive "Columns ..." blocks.
blocks = {};
block = [];
scale = 1;

for i = 1:numel(lines)
    line = lines{i};
    if startsWith(line, 'Columns ')
        if ~isempty(block)
            blocks{end + 1} = block; %#ok<AGROW>
            block = [];
        end
        continue
    end

    scale_match = regexp(line, ...
        '^([+-]?(?:\d+\.?\d*|\.\d+)[eE][+-]?\d+)\s*\*\s*$', ...
        'tokens', 'once');
    if ~isempty(scale_match)
        scale = str2double(scale_match{1});
        continue
    end

    row = sscanf(line, '%f').';
    if ~isempty(row)
        assert(isempty(block) || size(block, 2) == numel(row), ...
            'Inconsistent numeric row width.');
        block(end + 1, :) = row; %#ok<AGROW>
    end
end
if ~isempty(block)
    blocks{end + 1} = block;
end

assert(~isempty(blocks), 'No numeric data found.');
block_heights = cellfun(@(x) size(x, 1), blocks);
assert(numel(unique(block_heights)) == 1, ...
    'Wrapped matrix blocks have inconsistent heights.');
value = [blocks{:}];
value = scale * value;

% Logged vectors are printed vertically; preserve that column orientation.
if size(value, 1) == 1 && size(value, 2) > 1
    value = value.';
end
end
