Qw = 20;
Qp = 0;

Le = 0.1; % factor error between true and estimated L_delta

M = tf(216^2, [1, 158, 216^2]) * tf(1, [0.0126, 1]);

D = tf(1, 1, "iodelay", 0.01);

A = [0, 1; 0, 0];
B = [0; 1];
C = [1, 0];
Kp = sqrt(Qp);
Kw = sqrt(2*sqrt(Qp)+Qw);
K = -[Kp, Kw]*Le;
PC = ss([0,1;0,0]+B*K, B, C, 0);

G = M*PC*D;

margin(G)