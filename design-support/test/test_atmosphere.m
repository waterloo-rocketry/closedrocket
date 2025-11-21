clear

alt = 0:0.1:40; % altitude in km

for i = 1:length(alt)
    air_i = model_airdata(alt(i)*1000);
    p(i) = air_i.pressure;
    t(i) = air_i.temperature;
    rho(i) = air_i.density;
    mach(i) = air_i.mach;
end

%% Plot
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')

f_atm = figure(1);

tiledlayout(1,3,'TileSpacing','Compact');

plot_p = nexttile;
plot(p/1000, alt)
yline([max(alt), 32, 20, 11],':', {"Stratosphere 2", "Stratosphere 1", "Tropopause", "Troposphere"},'LabelVerticalAlignment','bottom', 'Interpreter','latex')
ylim([min(alt), max(alt)])
ylabel("Altitude $\ell, \: [\mathrm{km}]$")
xlabel("Pressure $p, \: [\mathrm{kPa}]$")
ax = gca;
ax.XGrid = 'on';
xticks(linspace(0, 100, 5));
xtickangle(0)
fontsize(14,"points")
 
plot_t = nexttile;
plot(t, alt)
yline([max(alt), 32, 20, 11],':', {"Strato. 2", "Stratosphere 1", "Tropopause", "Troposphere"},'LabelVerticalAlignment','bottom', 'Interpreter','latex')
ylim([min(alt), max(alt)])
ylabel("Altitude $\ell, \: [\mathrm{km}]$")
xlabel("Temperature $T, \: [\mathrm{K}]$")
ax = gca;
ax.XGrid = 'on';
xticks(linspace(200, 300, 3));
xtickangle(0)
fontsize(14,"points")

plot_rho = nexttile;
plot(rho, alt)
yline([max(alt), 32, 20, 11],':', {"Stratosphere 2", "Stratosphere 1", "Tropopause", "Troposphere"},'LabelVerticalAlignment','bottom', 'Interpreter','latex')
ylim([min(alt), max(alt)])
ylabel("Altitude $\ell, \: [\mathrm{km}]$")
xlabel("Density $\rho, \: [\mathrm{kg}/\mathrm{m}^3]$")
ax = gca;
ax.XGrid = 'on';
xticks(linspace(0, 1.2, 5));
xtickangle(0)
xlim([0, 1.25])
fontsize(14,"points")

% figure(4)
% plot(mach, alt)
% yline([max(alt), 32, 20, 11],'--', {"Stratosphere2", "Stratosphere", "Tropopause", "Troposphere"},'LabelVerticalAlignment','bottom', 'Interpreter','latex')
% ylim([min(alt), max(alt)])
% ylabel("Altitude $$l$$  [km]")
% xlabel("Speed of sound $$c$$ [m/s]")

%%

fontsize(f_atm, 12, "points")

exportgraphics(f_atm, 'plots/plot_atmosphere.pdf', 'ContentType', 'vector')

set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')

