function [u, r, C_l_delta] = controller_module(time, xR, pdyn, delta)
    % Top-level controller module.
    % u : (rad) control command, desired canard angle
    % r : (rad) roll angle target
    % time : (s) current time stamp
    % xR : [(rad) roll angle, (rad/s) roll rate] roll state
    % pdyn : (Pa) dynamic pressure

    persistent time_prev ctrl_mem

    if isempty(time_prev)
        time_prev = -0.01;
    end
    if isempty(ctrl_mem)
        ctrl_mem.coeffs = [2; 0];
        ctrl_mem.w = xR(2);
        ctrl_mem.P = diag([1e-5, 1e-9]);
        ctrl_mem.c_delta = 0;
        ctrl_mem.w_dot = 0;
    end

    time_launch = 20; % (s) pad delay time
    time_actuate = 10; % (s) time to actuation

    flight_time = time - time_launch;
    dt_ctrl = time - time_prev;

    [u, r, ctrl_mem] = controller_codegen_entry(flight_time, dt_ctrl, xR, pdyn, delta, ctrl_mem);

    if flight_time < time_actuate
        u = 0;
    end

    C_l_delta = ctrl_mem.coeffs(1);

    time_prev = time;
end
