function vec_norm = vec_norm(vec)
    if abs(norm(vec)) < 1e-6
        vec_norm = vec;
    else
        vec_norm = vec / norm(vec);
    end
end