q0 = [0.9 0 -0.1 0]';
q0 = q0/norm(q0);

w = 10*[10 2 5]';

dt = 0.05;

t_end = 10;

opts = odeset('RelTol',1e-8,'AbsTol',1e-12);
[t_ct, q_ct] = ode45(@(t,q) quat_ode(t,q,w), [0, t_end], q0, opts);

q_up = zeros(4, t_end/dt+1); 
q_up(:,1) = q0;
q_in = zeros(4, t_end/dt+1); 
q_in(:,1) = q0;

for t = 1:(t_end/dt)
    q_up(:,t+1) = quaternion_update(q_up(:,t), w, dt);
    q_in(:,t+1) = quaternion_increment(q_in(:,t),w,dt);
end

function qdot = quat_ode(t, q, w)
    qdot  = quaternion_derivative(q,w);
end

figure(1)
plot(t_ct, q_ct(:,1), DisplayName='ct')
hold on
plot(0:dt:t_end, q_up(1,:), DisplayName='up')
plot(0:dt:t_end, q_in(1,:), DisplayName='in')
legend("show")
hold off
