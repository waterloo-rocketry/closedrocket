Q_int = 0.2; Q_phi = 10; Q_omega = 10; R = 1;
Q = diag([Q_int Q_phi Q_omega]);

A = [0  -1   0;
     0   0   1;
     0   0   0];

B0 = [0;0;1];         % input is alpha directly
C  = [0 1 0];
E  = [1;0;0];

K = lqr(A,B0,Q,R)