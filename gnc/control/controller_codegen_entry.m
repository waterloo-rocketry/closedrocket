function [u, r, ctrl_mem] = controller_codegen_entry(time, dt_ctrl, xR, pdyn, delta, ctrl_mem)
    %#codegen

    % u : (rad) control command, desired canard angle
    % r : (rad) roll angle target
    % time : (s) current time stamp
    % xR : [(rad) roll angle, (rad/s) roll rate] reduced roll state
    % pdyn : (Pa) dynamic pressure

    persistent param
    if isempty(param)
        param = coder.load('gnc/model_params.mat');
    end

    %% Constants
    time_program = 15; % (s) time from launch to start of roll program
    u_max = deg2rad(20); % (rad) limit output to this angle
    L_min = 10; % (rad/s^2 / rad) limit roll control derivative for low authority conditions
    pdyn_min = 500; % (Pa) deactivate at low authority near apogee

    %% Reference signal
    %%% Generates reference signal r for roll program
    %%% includes multiple roll angle steps
    roll_step_times = [7, 12, 17, 24] + time_program;
    roll_step_targets = [0.5, -0.5, 0.5];

    r = 0;
    for step_idx = 1:(numel(roll_step_times) - 1)
        if time >= roll_step_times(step_idx) && time < roll_step_times(step_idx + 1)
            r = roll_step_targets(step_idx);
            break;
        end
    end
    % r = 0; % deactivate roll program

    %% controller algorithm
    %%% Computes control signal of the adaptive LQR controller.

    %%% Coefficient Estimation
    pdyn_params = pdyn * param.c_canard;

    ctrl_mem = controller_estimator(dt_ctrl, xR(2), delta, pdyn_params, ctrl_mem);

    C_l_delta = ctrl_mem.coeffs(1);
    L_delta = C_l_delta * pdyn_params / 2;

    if abs(L_delta) < L_min
        if L_delta >= 0
            L_delta = L_min;
        else
            L_delta = -L_min;
        end
    end

    %%% Control Law, Feedback + Feedforward tracking
    u = controller_law(xR, r, L_delta);

    %%% limit output to allowable angle
    u = min(max(u, -u_max), u_max); % upper bounds

    if pdyn < pdyn_min % disable during low control authority
        u = 0;
    end

end
