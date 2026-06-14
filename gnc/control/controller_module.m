function [cmd, ref, C_l_delta] = controller_module(time, xR, pdyn, encoder)
    % Top-level controller module.
    % cmd : (rad) control command, desired motor angle
    % ref : (rad) roll angle target
    % time : (s) current time stamp
    % xR : [(rad) roll angle, (rad/s) roll rate] roll state
    % pdyn : (Pa) dynamic pressure
    % encoder : (rad) motor angle measurement from encoder

    persistent time_prev ctrl_mem

    if isempty(time_prev)
        time_prev = -0.01;
    end
    if isempty(ctrl_mem)
        ctrl_mem.coeffs = [2; 0];
        ctrl_mem.w = xR(2);
        ctrl_mem.P = diag([1e-5, 1e-9]);
        ctrl_mem.delta_lp = 0;
        ctrl_mem.w_dot_lp = 0;
    end

    time_launch = 20; % (s) pad delay time
    time_actuate = 10; % (s) time to actuation

    flight_time = time - time_launch;
    dt_ctrl = time - time_prev;

    [cmd, ref, ctrl_mem] = controller_codegen_entry(flight_time, dt_ctrl, xR, pdyn, encoder, ctrl_mem);

    if flight_time < time_actuate
        cmd = 0;
    end

    C_l_delta = ctrl_mem.coeffs(1);

    time_prev = time;
end
