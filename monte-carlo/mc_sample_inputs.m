function sample = mc_sample_inputs(opts, baseline, run_index)
%MC_SAMPLE_INPUTS Sample one Monte Carlo input set.

    vars = struct();
    vars.engine_thrust_factor = clipped_normal( ...
        opts.performance.engine_thrust_factor);
    vars.act_backlash = clipped_normal( ...
        opts.performance.act_backlash_deg);

    canard_cant_zero_deg = clipped_normal( ...
        opts.geometry.canard_cant_zero_deg);
    vars.canard_cant_zero = deg2rad(canard_cant_zero_deg);

    canard_spine_effectiveness = clipped_normal( ...
        opts.geometry.canard_spine_effectiveness);
    vars.canard_spine_effectiveness = canard_spine_effectiveness;

    cp_shift_m = clipped_normal(opts.geometry.cp_shift_m);
    vars.nosecone_pos_x_cp = baseline.nosecone_pos_x_cp + cp_shift_m;
    vars.body_pos_x_cp = baseline.body_pos_x_cp + cp_shift_m;
    vars.tail_pos_x_cp = baseline.tail_pos_x_cp + cp_shift_m;
    vars.fin_pos_x_cp = baseline.fin_pos_x_cp + cp_shift_m;
    vars.fin_pos_x_cp_mach2 = baseline.fin_pos_x_cp_mach2 + cp_shift_m;

    radial_offset_m = clipped_normal( ...
        opts.geometry.canard_radial_cp_offset_m);
    vars.canard_radial_offset = ...
        baseline.canard_radial_offset + radial_offset_m;
    [vars.canard_pos_x, vars.canard_pos_r_mean, vars.canard_area] = ...
        aerosurface_canard( ...
        baseline.canard_chord_root, baseline.canard_height, ...
        baseline.canard_sweep_angle, baseline.canard_spine_width, ...
        canard_spine_effectiveness, vars.canard_radial_offset, ...
        baseline.canard_pos_x_roottip, baseline.rocket_diameter);
    vars.canard_pos_x = vars.canard_pos_x + cp_shift_m;

    reversal_occurs = rand() < opts.roll_reversal.probability;
    if reversal_occurs
        reversal_time_after_takeoff_s = uniform_sample( ...
            opts.roll_reversal.time_after_takeoff_s);
        reversal_factor = uniform_sample( ...
            opts.roll_reversal.factor_range);
        reversal_time_s = baseline.time_idle + ...
            reversal_time_after_takeoff_s;
        model_reversal_time_s = reversal_time_s;
    else
        reversal_time_after_takeoff_s = NaN;
        reversal_time_s = NaN;
        model_reversal_time_s = opts.stop_time + 1;
        reversal_factor = 1;
    end
    vars.roll_rev_grid_active = false;
    vars.canard_roll_reversal_time = model_reversal_time_s;
    vars.canard_roll_reversal_factor = reversal_factor;

    [wind_vars, wind_mode] = sample_wind(opts.wind);
    wind_names = fieldnames(wind_vars);
    for i = 1:numel(wind_names)
        name = wind_names{i};
        vars.(name) = wind_vars.(name);
    end

    sample.vars = vars;
    sample.meta.run_index = run_index;
    sample.meta.seed = opts.seed + run_index - 1;
    sample.meta.wind_mode = wind_mode;
    sample.parameters = struct( ...
        "run_index", run_index, ...
        "seed", sample.meta.seed, ...
        "wind_mode", wind_mode, ...
        "engine_thrust_factor", vars.engine_thrust_factor, ...
        "canard_roll_reversal_occurred", reversal_occurs, ...
        "canard_roll_reversal_time_after_takeoff_s", ...
        reversal_time_after_takeoff_s, ...
        "canard_roll_reversal_time_s", reversal_time_s, ...
        "canard_roll_reversal_factor", reversal_factor, ...
        "act_backlash_deg", vars.act_backlash, ...
        "canard_cant_zero_deg", canard_cant_zero_deg, ...
        "canard_spine_effectiveness", canard_spine_effectiveness, ...
        "cp_shift_m", cp_shift_m, ...
        "canard_radial_cp_offset_m", radial_offset_m);
end

function [vars, mode] = sample_wind(opts)
    mode = string(opts.mode);
    if mode == "mixed"
        if rand() < opts.historical_probability
            mode = "historical";
        else
            mode = "discrete";
        end
    end

    vars = struct();
    switch mode
        case "none"
            vars.wind_disturbance_enable = 0;
            vars.wind_discrete_enable = 1;

        case "historical"
            vars.wind_disturbance_enable = 1;
            vars.wind_discrete_enable = 0;
            [vars.wind_heights, vars.wind_vectors] = ...
                wind_historic(opts.layer_threshold_mps);

        case "discrete"
            vars.wind_disturbance_enable = 1;
            vars.wind_discrete_enable = 1;
            vars.wind_const_direction = 2 * pi * rand();
            vars.wind_const_strength = clipped_normal( ...
                opts.const_strength_mps);
            vars.wind_gust1_amplitude = clipped_normal( ...
                opts.gust_amplitude_mps);
            vars.wind_gust2_amplitude = clipped_normal( ...
                opts.gust_amplitude_mps);
            vars.wind_gust3_amplitude = clipped_normal( ...
                opts.gust_amplitude_mps);
            vars.wind_gust1_distribution = axis_distribution( ...
                [0.01; 0.9; 0.5]);
            vars.wind_gust2_distribution = axis_distribution( ...
                [-0.01; -0.5; -0.8]);
            vars.wind_gust3_distribution = axis_distribution( ...
                [0.01; 1; 1]);
    end
end

function value = clipped_normal(spec)
    lower_bound = min(spec.range);
    upper_bound = max(spec.range);
    standard_sample = randn();

    if standard_sample < 0
        sigma = (spec.mean - lower_bound) / 2.5;
    else
        sigma = (upper_bound - spec.mean) / 2.5;
    end

    value = spec.mean + sigma * standard_sample;
    value = max(min(value, upper_bound), lower_bound);
end

function value = uniform_sample(range)
    lower_bound = min(range);
    upper_bound = max(range);
    value = lower_bound + (upper_bound - lower_bound) * rand();
end

function distribution = axis_distribution(scale)
    distribution = (2 * rand(3, 1) - 1) .* scale;
end
