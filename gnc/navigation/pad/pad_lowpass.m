function filtered = pad_lowpass(filtered, sensor, alpha)
    %%% lowpass filter function used in pad filter
    if sensor.status == 1 
        filtered = alpha * sensor.meas + (1 - alpha) * filtered; 
    end
end
