function [u, r, C_l_delta] = controller_module(time, xR, pdyn, delta)
    % Top-level controller module.
    % u : control command, desired canard angle (rad)
    % r : roll angle target (rad)
    % time : current time stamp (s)    
    % xR : roll state [roll angle (rad), roll rate (rad/s)]
    % pdyn : dynamic pressure (Pa)

    time_launch = 20; % [s], pad delay time
    time_actuate = 10; % [s], time to actuation

    t = time - time_launch;

    [u, r, C_l_delta] = controller_codegen_entry(t, xR, pdyn, delta);

    if t < time_actuate
        u = 0;
    end
end

