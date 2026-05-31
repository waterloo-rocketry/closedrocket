function [u, r, C_l_delta] = controller_module(time, xR, pdyn, delta)
    % Top-level controller module.
    % u : control command, desired canard angle (rad)
    % r : roll angle target (rad)
    % time : current time stamp (s)    
    % xR : roll state [roll angle (rad), roll rate (rad/s)]
    % pdyn : dynamic pressure (Pa)

    persistent time_old coeffs w_old P_minus d_old w_dot_old

    if isempty(time_old)
        time_old = -0.01;
    end
    if isempty(coeffs)
        coeffs = [2; 0];
    end
    if isempty(w_old)
        w_old = xR(2);
    end
    if isempty(P_minus)
        P_minus = diag([1e-5, 1e-9]);
    end
    if isempty(d_old)
        d_old = 0;
    end
    if isempty(w_dot_old)
        w_dot_old = 0;
    end

    time_launch = 20; % [s], pad delay time
    time_actuate = 10; % [s], time to actuation

    t = time - time_launch;
    dt_ctrl = time - time_old;

    [u, r, coeffs_ret, w_old_ret, P_minus_ret, d_old_ret, w_dot_old_ret] = controller_codegen_entry(t, dt_ctrl, xR, pdyn, delta, w_old, coeffs, P_minus, d_old, w_dot_old);

    if t < time_actuate
        u = 0;
    end

    C_l_delta = coeffs(1);

    coeffs = coeffs_ret;
    P_minus = P_minus_ret;
    d_old = d_old_ret;
    w_old = w_old_ret;
    w_dot_old = w_dot_old_ret;
    time_old = time;
end