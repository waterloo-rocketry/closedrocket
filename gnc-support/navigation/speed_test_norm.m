A = rand(11);
A = A + 20*eye(11);
P = A'*A


tic
norm(P, 'fro');
toc
tic 
norm(P, 1);
toc
tic 
norm(P, inf);
toc
tic 
norm(P);
toc 