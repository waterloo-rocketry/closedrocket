function [u, r, C_l_delta] = controller_module(time, xR, pdyn, delta)
    % Top-level controller module.
    % u : control command, desired canard angle (rad)
    % r : roll angle target (rad)
    % time : current time stamp (s)    
    % xR : roll state [roll angle (rad), roll rate (rad/s)]
    % pdyn : dynamic pressure (Pa)

    %% settings
    time_launch = 10; % pad delay time
    time_coast = 10; % time from launch to burnout
    time_program = 20; % time from launch to start of roll program
    u_max = deg2rad(10); % limit output to this angle
    L_min = 10; % limit roll control derivative for low authority conditions

    %% parameters
    persistent param
    if isempty(param)
        param = coder.load("gnc/model_params.mat");
    end

    %% Reference signal
    % Generates reference signal r (rad) for roll program
    % includes multiple roll angle steps
    
    t = time - time_launch;
    r = 0;
    if t > (time_program + 7)
        if t < (time_program + 12)
            r = 0.5;
        elseif t < (time_program + 17)
            r = -0.5;
        elseif t < (time_program + 24)
            r = 0.5;
        elseif t > (time_program + 31)
            r = 0;
        end
    end

    %% controller algorithm
    % Computes control signal of the adaptive LQR controller.

    %%% Coefficient Estimation
    pdyn_params = pdyn * param.c_canard;
    C_l_delta = controller_estimator(time, xR(2), delta, pdyn_params);
       
    L_delta = C_l_delta * pdyn_params;
    L_delta = 1 / (max(min(1/L_delta, 1/L_min), -1/L_min)); % lower bounds

    %%% Feedback law
    u = 0;
    u = controller_law(xR, r, L_delta);

    %%% limit output to allowable angle
    u = min(max(u, -u_max), u_max); % upper bounds

    if t < time_coast % disable during boost
        u = 0;
    end
    if pdyn < 1000 % disable during low control auhtority at apogee
        u = 0;
    end
end

