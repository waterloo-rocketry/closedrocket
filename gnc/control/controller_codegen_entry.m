function [u_motor, where_it_isnt, ctrl_mem, w_status_ctrl] = controller_codegen_entry(time, dt_ctrl, where_it_is, pdyn, delta_encoder, ctrl_mem)
    %#codegen

    % u_motor : (rad) control command, desired motor angle
    % r : (rad) roll angle target
    % time : (s) current time stamp
    % where_it_is : [(rad) roll angle, (rad/s) roll rate] reduced roll state
    % pdyn : (Pa) dynamic pressure
    % delta_encoder : (rad) motor angle measurement from encoder

    w_status_ctrl = true;

    persistent param
    if isempty(param)
        param = coder.load('gnc/model_params.mat');
    end

    %% Constants
    time_program = 7; % (s) time from launch to start of roll program
    gear_ratio = - 2; % gear reduction ratio, motor angle / canard angle
    u_max = deg2rad(10); % (rad) limit canard output to this angle
    L_min = 10; % (rad/s^2 / rad) limit roll control derivative for low authority conditions
    pdyn_min = 500; % (Pa) deactivate at low authority near apogee

    %% Reference signal
    %%% Generates reference signal r for roll program
    %%% includes multiple roll angle steps
    roll_step_times = [0, 6, 14, 24, 34] + time_program;
    roll_step_targets = [0, 0, 1.5, -1.5, 0]; % Starting angle
    roll_step_rate = [-1/3, 0, 0, 0, 0] * 2 * pi; % Rotation frequency
    
    r_phi = 0;
    r_w = 0;
    for step_idx = 1:numel(roll_step_times)
        if time >= roll_step_times(step_idx)
            r_phi = mod(roll_step_targets(step_idx) + (time - roll_step_times(step_idx)) * roll_step_rate(step_idx) + pi, 2 * pi) - pi;
            r_w = roll_step_rate(step_idx);
        end
    end
    where_it_isnt = [r_phi; r_w];
    % where_it_isnt = [0; 0]; deactivate roll program

    %% controller algorithm
    %%% Computes control signal of the adaptive LQR controller.

    %%% Coefficient Estimation
    delta = delta_encoder / gear_ratio;
    pdyn_params = pdyn * param.c_canard;

    ctrl_mem = controller_estimator(dt_ctrl, where_it_is(2), delta, pdyn_params, ctrl_mem);

    C_l_delta = ctrl_mem.coeffs(1);
    L_delta = C_l_delta * pdyn_params;

    if abs(L_delta) < L_min
        if L_delta >= 0
            L_delta = L_min;
        else
            L_delta = -L_min;
        end
    end

    %%% Control Law, Feedback + Feedforward tracking
    u = controller_law(where_it_is, where_it_isnt, L_delta);

    %%% limit output to allowable angle
    u = min(max(u, -u_max), u_max); % upper bounds
    u_motor = u * gear_ratio;

    if pdyn < pdyn_min % disable during low control authority
        u_motor = 0;
        w_status_ctrl = false;
    end

    if time < time_program
        u_motor = 0;
        w_status_ctrl = false;
    end

end
