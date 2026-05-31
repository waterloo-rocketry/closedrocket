function [u, r, coeffs_ret, w_old_ret, P_minus_ret, d_old_ret, w_dot_old_ret] = controller_codegen_entry(time, dt_ctrl, xR, pdyn, delta, w_old, coeffs, P_minus, d_old, w_dot_old)
    %#codegen

    % u : control command, desired canard angle (rad)
    % r : roll angle target (rad)
    % time : current time stamp (s)    
    % xR : roll state [roll angle (rad), roll rate (rad/s)]
    % pdyn : dynamic pressure (Pa)

    persistent param
    if isempty(param)
        param = coder.load('gnc/model_params.mat');
    end

    %% Constants
    time_program = 15;
    u_max = deg2rad(10); % [rad], limit output to this angle
    L_min = 10; % [rad/s^2 (angular accelaration) / rad (canard angle)] limit roll control derivative for low authority conditions
    pdyn_min = 500; % [Pa] deactivate at low authority near apogee
    
    %% Reference signal
    %%% Generates reference signal r (rad) for roll program
    %%% includes multiple roll angle steps
    if time >= (time_program + 7) && time < (time_program + 12)
        r = 0.5;
    elseif time >= (time_program + 12) && time < (time_program + 17)
        r = -0.5;
    elseif time >= (time_program + 17) && time < (time_program + 24)
        r = 0.5;
    else
        r = 0;
    end
    % r = 0; % deactivate roll program

    %% controller algorithm
    %%% Computes control signal of the adaptive LQR controller.
  
    %%% Coefficient Estimation
    pdyn_params = pdyn * param.c_canard;

    [coeffs_ret, w_old_ret, P_minus_ret, d_old_ret, w_dot_old_ret] = controller_estimator(dt_ctrl, xR(2), delta, pdyn_params, w_old, coeffs, P_minus, d_old, w_dot_old);
    
    C_l_delta = coeffs_ret(1);
    L_delta = C_l_delta * pdyn_params;
    
    if abs(L_delta) < L_min
        if L_delta >= 0
            L_delta = L_min;
        else
            L_delta = -L_min;
        end
    end
   
    %%% Control Law, Feedback + Feedforward tracking
    u = 0;
    u = controller_law(xR, r, L_delta);

    %%% limit output to allowable angle
    u = min(max(u, -u_max), u_max); % upper bounds

    if pdyn < pdyn_min % disable during low control authority
        u = 0;
    end

end