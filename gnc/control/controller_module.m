function [u, r, C_l_delta] = controller_module(time, xR, pdyn, delta)
    % Top-level controller module.
    % u : control command, desired canard angle (rad)
    % r : roll angle target (rad)
    % time : current time stamp (s)    
    % xR : roll state [roll angle (rad), roll rate (rad/s)]
    % pdyn : dynamic pressure (Pa)

    persistent time_old ctrl_mem_in

    if isempty(time_old)
        time_old = -0.01;
    end
    if isempty(ctrl_mem_in)
        ctrl_mem_in.coeffs = [2; 0];
        ctrl_mem_in.w_old = xR(2);
        ctrl_mem_in.P_minus = diag([1e-5, 1e-9]);
        ctrl_mem_in.d_old = 0;
        ctrl_mem_in.w_dot_old = 0;
    end
    % 
    % if isempty(coeffs)
    %     coeffs = [2; 0];
    % end
    % if isempty(w_old)
    %     w_old = xR(2);
    % end
    % if isempty(P_minus)
    %     P_minus = diag([1e-5, 1e-9]);
    % end
    % if isempty(d_old)
    %     d_old = 0;
    % end
    % if isempty(w_dot_old)
    %     w_dot_old = 0;
    % end

    time_launch = 20; % [s], pad delay time
    time_actuate = 10; % [s], time to actuation

    flight_time = time - time_launch;
    dt_ctrl = time - time_old;

    [u, r, ctrl_mem_out] = controller_codegen_entry(flight_time, dt_ctrl, xR, pdyn, delta, ctrl_mem_in);

    if flight_time < time_actuate
        u = 0;
    end

    C_l_delta = ctrl_mem_out.coeffs(1);

    ctrl_mem_in = ctrl_mem_out;
    time_old = time;
end