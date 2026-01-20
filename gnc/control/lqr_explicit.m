syms L R real
syms q1 q2 real

assume(L > 0)
assume(R > 0)
assume(q1 > 0)
assume(q1 > 0)

% R = 1
% q1 = 2; q2 = 3;

A = [0, 1; 0, 0]
B = [0; L]

Q = [q1, 0; 0, q2]

H = [A, -B/R*B'; -Q, -A']

[vec, val] = eig(H)
vec = simplify(expand(vec))
val = simplify(expand(val))

U = vec; % stable eigenmodes on left side

% U(:,2) = vec(:,3); U(:,3) = vec(:,2) % in case stable eigenvalues are 1, 3 (not 1, 2)
% U(:,2) = vec(:,4); U(:,4) = vec(:,2)
% U(:,1) = vec(:,3); U(:,3) = vec(:,1)

P = U(3:4,1:2) * inv(U(1:2,1:2))

K = simplify(expand(-1/R*B'*P))


L = 10;
q1 = 1;
q2 = 1;
R = 10;

double(subs(K))
double(subs(val))
double(subs(U))


