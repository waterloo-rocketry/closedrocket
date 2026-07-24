%% Logged testcase MATLAB versus SIL equivalence tests

script_dir = fileparts(mfilename('fullpath'));
testcase_file = fullfile(script_dir, 'testcases.mat');
assert(isfile(testcase_file), ...
    'Missing %s. Run generate_testcases_mat first.', testcase_file);
load(testcase_file, 'controller_testcases', 'navigator_testcases');

absTol = 1e-10;
relTol = 1e-9;

%% Controller testcases

for testcase_index = 1:numel(controller_testcases)
    testcase = controller_testcases(testcase_index);

    [u_m, r_m, ctrl_mem_m, w_status_ctrl_m] = ...
        controller_codegen_entry(testcase.time, testcase.dt, ...
        testcase.roll_state, testcase.pdyn, testcase.delta, ...
        testcase.ctrl_mem);

    [u_s, r_s, ctrl_mem_s, w_status_ctrl_s] = ...
        GNC_codegen_SIL_sil('controller_codegen_entry', ...
        testcase.time, testcase.dt, testcase.roll_state, testcase.pdyn, ...
        testcase.delta, testcase.ctrl_mem);

    prefix = sprintf('controller testcase %d', testcase_index);
    compare_value([prefix '.u'], u_m, u_s, absTol, relTol);
    compare_value([prefix '.r'], r_m, r_s, absTol, relTol);
    compare_value([prefix '.ctrl_mem'], ctrl_mem_m, ctrl_mem_s, absTol, relTol);
    compare_value([prefix '.w_status_ctrl'], w_status_ctrl_m, ...
        w_status_ctrl_s, absTol, relTol);
end

%% Navigation testcases

for testcase_index = 1:numel(navigator_testcases)
    testcase = navigator_testcases(testcase_index);

    [x_m, P_m, bias_m, sens_filt_m, cov_norm_m, roll_state_m, ...
        pdyn_nav_m, w_status_nav_m] = navigation_codegen_entry( ...
        testcase.dt, testcase.flight_phase, testcase.x, testcase.P, ...
        testcase.bias, testcase.sens_filt, testcase.sens_input);

    [x_s, P_s, bias_s, sens_filt_s, cov_norm_s, roll_state_s, ...
        pdyn_nav_s, w_status_nav_s] = GNC_codegen_SIL_sil( ...
        'navigation_codegen_entry', testcase.dt, testcase.flight_phase, ...
        testcase.x, testcase.P, testcase.bias, testcase.sens_filt, ...
        testcase.sens_input);

    prefix = sprintf('navigation testcase %d', testcase_index);
    compare_value([prefix '.x'], x_m, x_s, absTol, relTol);
    compare_value([prefix '.P'], P_m, P_s, absTol, relTol);
    compare_value([prefix '.bias'], bias_m, bias_s, absTol, relTol);
    compare_value([prefix '.sens_filt'], sens_filt_m, sens_filt_s, absTol, relTol);
    compare_value([prefix '.cov_norm'], cov_norm_m, cov_norm_s, absTol, relTol);
    compare_value([prefix '.roll_state'], roll_state_m, roll_state_s, absTol, relTol);
    compare_value([prefix '.pdyn'], pdyn_nav_m, pdyn_nav_s, absTol, relTol);
    compare_value([prefix '.w_status_nav'], w_status_nav_m, ...
        w_status_nav_s, absTol, relTol);
end

fprintf('\nAll %d controller and %d navigation MATLAB vs SIL comparisons passed.\n', ...
    numel(controller_testcases), numel(navigator_testcases));

clear GNC_codegen_SIL_sil;


%% Local comparison helper

function compare_value(name, a, b, absTol, relTol)
if isstruct(a)
    fields_a = fieldnames(a);
    fields_b = fieldnames(b);

    assert(isequal(sort(fields_a), sort(fields_b)), ...
        'Field mismatch in %s', name);

    for i = 1:numel(fields_a)
        field = fields_a{i};
        compare_value([name '.' field], a.(field), b.(field), absTol, relTol);
    end

elseif isnumeric(a)
    assert(isequal(size(a), size(b)), 'Size mismatch in %s', name);

    difference = abs(a - b);
    tolerance = absTol + relTol * max(abs(a), abs(b));

    if any(difference(:) > tolerance(:))
        error('Numeric mismatch in %s. Max error = %.3e', ...
            name, max(difference(:)));
    end

elseif islogical(a)
    assert(isequal(a, b), 'Logical mismatch in %s', name);

else
    assert(isequal(a, b), 'Value mismatch in %s', name);
end
end
