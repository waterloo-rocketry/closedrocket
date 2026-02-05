A = [0,1;0,0];
B  = [0; 1];

L = 100:100:1000;

Qp_t = 1;
Qw_t = 2;

Qp_a = 1;
Qw_a = 2;

for i = 1:length(L)
    [Kt, Ktr] = controller_gains_torque(L(i), Qp_t, Qw_t)
    [Ka, Kar] = controller_gains_angle(L(i), Qp_a/L(i)^2, Qw_a/L(i)^2)

    sys_t = ss(A, B*L(i), Kt, 0);
    sys_a = ss(A, B*L(i), Ka, 0);
    sys_array_t(:,:,1,i) = sys_t;
    sys_array_a(:,:,1,i) = sys_a;

    cl_t = Ktr*ss(A+B*L(i)*Kt, B*L(i), eye(2), 0);
    cl_a = Kar*ss(A+B*L(i)*Ka, B*L(i), eye(2), 0);
    cl_array_t(:,:,1,i) = cl_t;
    cl_array_a(:,:,1,i) = cl_a;

end

figure(1)
bode(sys_array_a, "r")
hold on
bode(sys_array_t, "b")
hold on

figure(2)
step(cl_array_a, "r")
hold on
step(cl_array_t, "b")
hold on