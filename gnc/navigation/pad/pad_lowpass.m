function filtered = pad_lowpass(filtered, measured, alpha)
    %%% lowpass filter function used in pad filter
    if isempty(filtered) % initialize
        filtered = measured;  
    else 
        filtered = alpha * measured + (1 - alpha) * filtered; 
    end
end