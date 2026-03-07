function [q_normed] = quaternion_norm(q_unnormed)
    % norms quaternion

    %%% inverse quaternion 
    q_mag = norm(q_unnormed);
    q_normed = q_unnormed / q_mag;
end