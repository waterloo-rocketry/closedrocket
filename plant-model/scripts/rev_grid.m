function roll_reversal_factor = rev_grid(v, a) 
    
    x = rad2deg(a);
    
    b = 1/((x-3)^2+1);
    c = 1/((x+3)^2+1);
    d = exp(-x^2/100);
    e = 10*(x/(10+log(v+1)));
    
    roll_reversal_factor = 3 * (b + c) * d * (sin(e)^2 + 2)/2 * cos(e);
end