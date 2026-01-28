 load("controller/gains.mat", "Ks", "P_mesh", "C_mesh");

% samplep = 1e5; samplec = 1.5;
% for i=1:4
%     K(i) = interp2(P_mesh, C_mesh, Ks(:,:,i), samplec, samplep, 'linear');
% end

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')

f_gains = figure(4);

subplot(2,2,1)
% [P_mesh,C_mesh] = meshgrid(Cls,Ps);
surfl(P_mesh,C_mesh,Ks(:,:,1), 'FaceAlpha',0.5, 'EdgeColor','none')
hold on
% scatter3(samplec, samplep ,K(1), 20, "k", "o", "filled")
hold off
xlabel("$C_{L_\delta}$")
ylabel("$\bar p \; [\mathrm{kPa}]$")
yticks = get(gca, 'YTick');
set(gca, 'YTickLabel', yticks / 1e3);
zlabel("$K_\phi$")
% zlim([-1,1])
ax = gca;
ax.XLabel.Position = [0 -1e5 -0.3];   
ax.YLabel.Position = [-30 3e5 -0.3]; 
ax.ZLabel.Position = [-30 6.5e5 0]; 

% figure(2)
subplot(2,2,2)
% [P_mesh,C_mesh] = meshgrid(Cls,Ps);
surfl(P_mesh,C_mesh,Ks(:,:,2), 'FaceAlpha',0.5, 'EdgeColor','none')
hold on
% scatter3(samplec, samplep ,K(2), 20, "k", "o", "filled")
hold off
xlabel("$C_{L_\delta}$")
ylabel("$\bar p \; [\mathrm{kPa}]$")
yticks = get(gca, 'YTick');
set(gca, 'YTickLabel', yticks / 1e3);
zlabel("$K_{\omega_x}$")
zlim([-1,1])
ax = gca;
ax.XLabel.Position = [0 -1e5 -1.5];   
ax.YLabel.Position = [-30 3e5 -1.5]; 
ax.ZLabel.Position = [-30 6.5e5 0]; 

% figure(3)
subplot(2,2,3)
% [P_mesh,C_mesh] = meshgrid(Cls,Ps);
surfl(P_mesh,C_mesh,Ks(:,:,3), 'FaceAlpha',0.5, 'EdgeColor','none')
hold on
% scatter3(samplec, samplep ,K(3), 20, "k", "o", "filled")
hold off
xlabel("$C_{L_\delta}$")
ylabel("$\bar p \; [\mathrm{kPa}]$")
yticks = get(gca, 'YTick');
set(gca, 'YTickLabel', yticks / 1e3);
zlabel("$K_\delta$")
% zlim([-4,0])ax = gca;
ax = gca;
ax.XLabel.Position = [0 -1e5 -5];   
ax.YLabel.Position = [-30 3e5 -5]; 
ax.ZLabel.Position = [-30 6.5e5 -2]; 

% figure(4)
subplot(2,2,4)
% [P_mesh,C_mesh] = meshgrid(Cls,Ps);
surfl(P_mesh,C_mesh,Ks(:,:,3), 'FaceAlpha',0.5, 'EdgeColor','none')
hold on
% scatter3(samplec, samplep ,K(4), 20, "k", "o", "filled")
hold off
xlabel("$C_{L_\delta}$")
ylabel("$\bar p \; [\mathrm{kPa}]$")
yticks = get(gca, 'YTick');
set(gca, 'YTickLabel', yticks / 1e3);
zlabel("$K_r$")
% zlim([-1,1])
ax = gca;
ax.XLabel.Position = [0 -1e5 -5];   
ax.YLabel.Position = [-30 3e5 -5]; 
ax.ZLabel.Position = [-30 6.5e5 -2]; 

fontsize(f_gains, 12, "points")

exportgraphics(f_gains, 'plots/plot_gains.pdf', 'ContentType','vector' )

set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')
