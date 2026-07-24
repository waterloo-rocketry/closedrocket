function vec_norm = vec_norm(vec)
    if norm(vec) == 0
        vec_norm = vec;
    else
        vec_norm = vec / norm(vec);
    end
end