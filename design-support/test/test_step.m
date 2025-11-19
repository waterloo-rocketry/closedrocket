clear
%% initials
CL = 4; % canard coefficient
alt = 1000; % altitude for dyn pressure
rho = model_airdata(alt).density;

V = linspace(30, 900, 20);
P = 0.5 * rho * V.^2;
T_sample = 0.005; % sampling time of the loop

clear control_scheduler
clear model_roll

%% test lqr + step
for i=1:length(P)   
    Ks = control_scheduler([P(i), CL]);
    K_pre = Ks(4);
    K = Ks(1:3);

    [A, B, C, ~] = model_roll(P(i), CL);
    sys_plant = c2d(ss(A, B, eye(3), 0), T_sample, 'zoh');
    [phi, gamma] = ssdata(sys_plant);

    %%% rolloff filter
    f_rolloff = 100; % [rad/s] rolloff frequency
    lowpass = c2d(tf(f_rolloff, [1, f_rolloff]), T_sample);

    sys_ol = - K_pre * ss(phi, gamma, K, 0, T_sample);
    % sys_ol = - K_pre * ss(A, B, K, 0);
    sys_cl = K_pre * ss(phi+gamma*K, gamma, eye(3), 0, T_sample);
    % sys_cl = K_pre * ss(A+B*K, B, eye(3), 0);
    

    sys_array(:,:,1,i) = sys_cl;
    sys_array_open(:,:,1,i) = sys_ol;

    if i == 1
        sys_min = sys_cl;
    elseif i == length(P)
        sys_max = sys_cl;
    end
    % step(sys_cl, steptime);
    % hold on
end

%% Figure
load("plots\mumColors.mat")

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')


t = 0:T_sample:3;

u = ones(length(t),1);
yout = zeros(length(t), 3, length(P));
for i = 1:length(P)
    yout(:,:,i) = lsim(sys_array(:,:,1,i), u, t);
end

f_step = figure(1);
names = ["$\phi$", "$\omega_x$", "$\delta$"];
for k = 1:3
    subplot(3,1,k)
    for i = 2:length(P)-1
        plot(t, yout(:,k,i), "Color", col.blue);
        hold on
    end
    plot(t, yout(:,k,1), "Color", col.green, "LineWidth",1);
    plot(t, yout(:,k,length(P)), "Color", col.red, "LineWidth",1);
    xlabel("Time [s]")
    ylabel(names(k))
    grid on
    hold off
end
fontsize(f_step, 12, "points")


w = logspace(-2,3,1000)';
for i = 1:length(P)
    [gain{i}, phase{i}] = bode(sys_array(1,1,1,i),w);
end

f_bode = figure(2);
subplot(2,1,1)
for i = 2:length(P)-1
    g = gain{i};
    g = squeeze(g(1,1,:));
    semilogx(w, 20*log10(g), "Color", col.blue)
    hold on   
end
g = gain{1}; g = squeeze(g(1,1,:));
semilogx(w, 20*log10(g), "Color", col.green, "LineWidth",1)
g = gain{end}; g = squeeze(g(1,1,:));
semilogx(w, 20*log10(g), "Color", col.red, "LineWidth",1)
xlim([1e-2, pi/T_sample])
xlabel("Frequency [rad/s]")
ylabel("Gain [dB]")
ylim([-100, 0])
yticks(-100:20:0)
grid on
hold off
subplot(2,1,2)
for i = 1:length(P)
    p = phase{i};
    p = squeeze(p(1,1,:));
    semilogx(w, p, "Color", col.blue)
    hold on   
end
p = phase{1}; p = squeeze(p(1,1,:));
semilogx(w, p, "Color", col.green, "LineWidth",1)
p = phase{end}; p = squeeze(p(1,1,:));
semilogx(w, p, "Color", col.red, "LineWidth",1)
xlim([1e-2, pi/T_sample])
ylim([-360, 0])
yticks(-360:90:0)
xlabel("Frequency [rad/s]")
ylabel("Phase [deg]")
grid on
hold off
fontsize(f_bode, 12, "points")


w = logspace(-2,3,1000)';
for i = 1:length(P)
    [gain{i}, phase{i}] = bode(sys_array_open(1,1,1,i),w);
end

f_bode_open = figure(3);
subplot(2,1,1)
for i = 2:length(P)-1
    g = gain{i};
    g = squeeze(g(1,1,:));
    semilogx(w, 20*log10(g), "Color", col.blue)
    hold on   
end
g = gain{1}; g = squeeze(g(1,1,:));
semilogx(w, 20*log10(g), "Color", col.green, "LineWidth",1)
g = gain{end}; g = squeeze(g(1,1,:));
semilogx(w, 20*log10(g), "Color", col.red, "LineWidth",1)
yline(0, "k-")
xlim([1e-2, pi/T_sample])
xlabel("Frequency [rad/s]")
ylabel("Gain [dB]")
ylim([-60, 60])
yticks(-60:20:60)
grid on
hold off

subplot(2,1,2)
for i = 1:length(P)
    p = phase{i};
    p = squeeze(p(1,1,:));
    semilogx(w, p, "Color", col.blue)
    hold on   
end
p = phase{1}; p = squeeze(p(1,1,:));
semilogx(w, p, "Color", col.green, "LineWidth",1)
p = phase{end}; p = squeeze(p(1,1,:));
semilogx(w, p, "Color", col.red, "LineWidth",1)
xlim([1e-2, pi/T_sample])
ylim([-180, -90])
yticks(-180:45:-90)
xlabel("Frequency [rad/s]")
ylabel("Phase [deg]")
grid on
hold off
fontsize(f_bode_open, 12, "points")

exportgraphics(f_step, 'plots/plot_step.pdf', 'ContentType', 'vector')
exportgraphics(f_bode, 'plots/plot_bode.pdf', 'ContentType', 'vector')
exportgraphics(f_bode_open, 'plots/plot_bode_open.pdf', 'ContentType', 'vector')

set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')
