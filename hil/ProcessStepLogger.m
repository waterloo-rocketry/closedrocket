load wall_step_log.mat

figure;
plot(sim_time, ratio);
grid on;
xlabel("Simulation time [s]");
ylabel("Wall dt / expected Ts");
title("Per-step real-time factor");
yline(1, "--");

slow_idx = ratio > 1.0;

Tslow = table( ...
    sim_time(slow_idx), ...
    wall_dt(slow_idx), ...
    ratio(slow_idx), ...
    'VariableNames', ["sim_time_s", "wall_dt_s", "real_time_ratio"]);

Tslow = sortrows(Tslow, "real_time_ratio", "descend");

disp(Tslow(1:min(30,height(Tslow)), :));
clear sim_time wall_dt ratio Ts;
run("configure_plant_model");