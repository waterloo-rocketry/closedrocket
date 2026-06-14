L = rand(3); L = L*L';
H = rand(3,11);
P = rand(11)*5; P = P*P';

tic
K1 = P * H' *inv(L);
toc

tic 
K2 = P * H' / L;
toc