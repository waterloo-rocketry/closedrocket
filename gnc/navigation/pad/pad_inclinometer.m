function [q] = pad_inclinometer(a)
    %%% computes initial orientation of stationary body from gravity acceleration
    %%% Input: 3D acceleration vector
    %%% Output: Orientation quaternion

    %%% normed gravity vector in body-fixed frame
    a_norm = norm(a);
    if a_norm < 1e-6
        q = [1; 0; 0; 0]; % fallback: assume exactly vertical
        return;
    end
    A = a / a_norm; % unit vector of gravity direction

    %%% determine initial orientation quaternion
    qw = sqrt( 0.5 + 0.5*A(1) );
    qx = 0;
    if qw == 0 % exact upside down case
        qy = 1; % either qy = 1 or qz = 1, this is arbitrary
        qz = 0;
    else
        qy = 0.5 * A(3) / qw;
        qz = -0.5 * A(2) / qw;
    end
    q = [qw; qx; qy; qz];
    q = q / norm(q);
end
