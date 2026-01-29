function [params, K] = param_estimator_gradient(time, w, d)
    % w : angular rate measurement
    % d : canard angle measurement / command
    % p : parameters
    % r : regressors
    % Q : weights
    % K : gain

    
    %%% initialize persistent variables
    persistent t p w_old
    if isempty(t)
        t = -0.01;
    end
    if isempty(p)
        p = [-0.5; 0.2]; % parameters, initial guess of CL_delta and CL_0
    end
    if isempty(w_old)
        w_old = w;
    end 
    
    %%% set up 
    dt = time - t;
    dw = (w - w_old) / dt; % actual roll acceleration  
    r = [d; 1]; % regressors, delta (cmd or encoder) and constant disturbance

    %%% measurement prediction
    pdyn_area_arm_inertia = 10; % everything in L except coefficients. match with gain in sim
    dw_hat = r' * p * pdyn_area_arm_inertia;

    %%% correction
    Q = diag([ 2, 1]) * 0.001; % weights
    Q = Q / (norm(r)^2 + 0.1); % normalized gradient
    K = Q * r; % gain
    p = p + K * (dw - dw_hat); % correct parameters
    
    %%% outputs and update persistent
    t = time;
    w_old = w;
    params = p;
end
