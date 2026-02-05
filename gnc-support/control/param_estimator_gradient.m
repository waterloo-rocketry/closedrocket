function [params, K] = param_estimator_gradient(time, w, d, c)
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
        p = [2; 0]; % parameters, initial guess of CL_delta and CL_0
    end
    if isempty(w_old)
        w_old = w;
    end 
    
    %%% set up 
    dt = time - t;
    dw = (w - w_old) / dt; % actual roll acceleration  
    r = [d; 0.1] ; % regressors, delta (cmd or encoder) and constant disturbance
    r = r * c;

    %%% measurement prediction
    dw_hat = r' * p;

    %%% correction
    Q = diag([1, 0.1]); % weights
    Q = Q / (norm(r)^2); % normalized gradient
    K = Q * r; % gain
    p = p + K * (dw - dw_hat); % correct parameters
    
    %%% outputs and update persistent
    t = time;
    w_old = w;
    params = p;
end
