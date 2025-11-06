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

f_p = figure(1);
plot(p/1000, alt)
yline([max(alt), 32, 20, 11],'--', {"Stratosphere 2", "Stratosphere 1", "Tropopause", "Troposphere"},'LabelVerticalAlignment','bottom', 'Interpreter','latex')
ylim([min(alt), max(alt)])
ylabel("Altitude $\ell, \: [\mathrm{km}]$")
xlabel("Pressure $p, \: [\mathrm{kPa}]$")
ax = gca;
ax.XGrid = 'on';
xticks(0:25:100);
fontsize(14,"points")
 
f_t = figure(2);
plot(t, alt)
yline([max(alt), 32, 20, 11],'--', {"Stratosphere 2", "Stratosphere 1", "Tropopause", "Troposphere"},'LabelVerticalAlignment','bottom', 'Interpreter','latex')
ylim([min(alt), max(alt)])
ylabel("Altitude $\ell, \: [\mathrm{km}]$")
xlabel("Temperature $T, \: [\mathrm{K}]$")
ax = gca;
ax.XGrid = 'on';
xticks(200:25:300);
fontsize(14,"points")

f_d = figure(3);
plot(rho, alt)
yline([max(alt), 32, 20, 11],'--', {"Stratosphere 2", "Stratosphere 1", "Tropopause", "Troposphere"},'LabelVerticalAlignment','bottom', 'Interpreter','latex')
ylim([min(alt), max(alt)])
ylabel("Altitude $\ell, \: [\mathrm{km}]$")
xlabel("Density $\rho, \: [\mathrm{kg}/\mathrm{m}^3]$")
ax = gca;
ax.XGrid = 'on';
xticks(0:0.4:1.5);
fontsize(14,"points")

% figure(4)
% plot(mach, alt)
% yline([max(alt), 32, 20, 11],'--', {"Stratosphere2", "Stratosphere", "Tropopause", "Troposphere"},'LabelVerticalAlignment','bottom', 'Interpreter','latex')
% ylim([min(alt), max(alt)])
% ylabel("Altitude $$l$$  [km]")
% xlabel("Speed of sound $$c$$ [m/s]")

%%

fontsize(f_p, 12, "points")
fontsize(f_t, 12, "points")
fontsize(f_d, 12, "points")

exportgraphics(f_p, 'plots/plot_atmosphere-pressure.pdf')
exportgraphics(f_t, 'plots/plot_atmosphere-temperature.pdf')
exportgraphics(f_d, 'plots/plot_atmosphere-density.pdf')

set(groot, 'defaultAxesTickLabelInterpreter','remove')
set(groot, 'defaultLegendInterpreter','remove')
set(groot, 'DefaultTextInterpreter', 'remove')

