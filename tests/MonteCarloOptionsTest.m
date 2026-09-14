classdef MonteCarloOptionsTest < matlab.unittest.TestCase
    %MONTECARLOOPTIONSTEST Fast tests for Monte Carlo configuration/sampling.

    properties
        ProjectRoot
    end

    methods (TestClassSetup)
        function addProjectPaths(testCase)
            testCase.ProjectRoot = fileparts( ...
                fileparts(mfilename("fullpath")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testCase.ProjectRoot, "monte-carlo")));
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture( ...
                fullfile(testCase.ProjectRoot, "plant-model", "scripts")));
        end
    end

    methods (TestMethodSetup)
        function seedRandomStream(testCase)
            original_rng = rng;
            testCase.addTeardown(@() rng(original_rng));
            rng(42, "twister");
        end
    end

    methods (Test, TestTags = {'Unit'})
        function testDefaultsMatchAgreedParameterRanges(testCase)
            opts = mc_default_options();
            geometry_fields = string(fieldnames(opts.geometry));

            testCase.verifyTrue(ismember( ...
                "canard_spine_effectiveness", geometry_fields));
            testCase.verifyTrue(ismember("cp_shift_m", geometry_fields));
            testCase.verifyTrue(ismember( ...
                "canard_radial_cp_offset_m", geometry_fields));
            testCase.verifyFalse(ismember( ...
                "canard_pitch_gain", geometry_fields));
            testCase.verifyFalse(ismember( ...
                "canard_pitch_zero_deg", geometry_fields));
            testCase.verifyFalse(ismember( ...
                "cp_canard_shift_m", geometry_fields));
            testCase.verifyFalse(ismember( ...
                "cp_common_shift_m", geometry_fields));
            testCase.verifyFalse(ismember( ...
                "cp_fin_shift_m", geometry_fields));
            testCase.verifyFalse(ismember( ...
                "fin_cant_total_deg", geometry_fields));
            testCase.verifyEqual( ...
                opts.performance.engine_thrust_factor.range, [0.8, 1.1]);
            testCase.verifyEqual( ...
                opts.performance.act_backlash_deg.range, [0, 0.5]);
            testCase.verifyEqual( ...
                opts.geometry.canard_cant_zero_deg.range, [0, 2]);
            testCase.verifyEqual( ...
                opts.geometry.canard_spine_effectiveness.range, [0.3, 1]);
            testCase.verifyEqual( ...
                opts.geometry.cp_shift_m.range, [-0.1, 0.1]);
            testCase.verifyEqual( ...
                opts.geometry.canard_radial_cp_offset_m.range, [-0.01, 0.02]);
            testCase.verifyEqual(opts.roll_reversal.probability, 0.5);
            testCase.verifyEqual(opts.roll_reversal.factor_range, [-1, 2]);
            testCase.verifyEqual( ...
                opts.roll_reversal.time_after_takeoff_s, [0, 40]);
        end

        function testSpineEffectivenessRecomputesCanardGeometry(testCase)
            opts = fixed_options();
            baseline = geometry_baseline();

            sample = mc_sample_inputs(opts, baseline, 1);

            testCase.verifyEqual( ...
                sample.vars.canard_area, 0.0324, AbsTol=1e-12);
            testCase.verifyEqual( ...
                sample.vars.canard_pos_x, -1.047407407407407, ...
                AbsTol=1e-12);
            testCase.verifyEqual( ...
                sample.vars.canard_pos_r_mean, 0.185987654320988, ...
                AbsTol=1e-12);
            testCase.verifyEqual( ...
                sample.vars.nosecone_pos_x_cp, 0.11, AbsTol=1e-12);
            testCase.verifyEqual( ...
                sample.vars.fin_pos_x_cp, -3.99, AbsTol=1e-12);
            testCase.verifyEqual( ...
                sample.vars.canard_radial_offset, 0.025, AbsTol=1e-12);
            testCase.verifyEqual( ...
                rad2deg(sample.vars.canard_cant_zero), 0.3, AbsTol=1e-12);
            testCase.verifyEqual( ...
                sample.parameters.canard_spine_effectiveness, 0.4);
            testCase.verifyTrue( ...
                sample.parameters.canard_roll_reversal_occurred);
            testCase.verifyEqual( ...
                sample.vars.canard_roll_reversal_time, 38);
            testCase.verifyEqual( ...
                sample.vars.canard_roll_reversal_factor, -0.5);
            testCase.verifyFalse(isfield(sample.vars, "act_gear_ratio"));
        end

        function testSamplingIsReproducibleForASeed(testCase)
            opts = mc_default_options();
            opts.wind.mode = "none";
            baseline = geometry_baseline();

            rng(1234, "twister");
            first_sample = mc_sample_inputs(opts, baseline, 7);
            rng(1234, "twister");
            second_sample = mc_sample_inputs(opts, baseline, 7);

            testCase.verifyEqual(second_sample, first_sample);
        end

        function testNoReversalUsesFiniteModelSwitchTime(testCase)
            opts = fixed_options();
            opts.roll_reversal.probability = 0;
            baseline = geometry_baseline();

            sample = mc_sample_inputs(opts, baseline, 1);

            testCase.verifyFalse( ...
                sample.parameters.canard_roll_reversal_occurred);
            testCase.verifyTrue(isnan( ...
                sample.parameters.canard_roll_reversal_time_s));
            testCase.verifyEqual( ...
                sample.vars.canard_roll_reversal_time, opts.stop_time + 1);
            testCase.verifyEqual( ...
                sample.vars.canard_roll_reversal_factor, 1);
        end

        function testInvalidWindModeIsRejected(testCase)
            opts = mc_default_options();
            opts.wind.mode = "tornado";

            testCase.verifyError(@() mc_run_sweep(opts), ...
                "mc_run_sweep:InvalidWindMode");
        end
    end
end

function opts = fixed_options()
    opts = mc_default_options();
    opts.wind.mode = "none";
    opts.geometry.canard_cant_zero_deg = fixed_spec(0.3);
    opts.geometry.canard_spine_effectiveness = fixed_spec(0.4);
    opts.geometry.cp_shift_m = fixed_spec(0.01);
    opts.geometry.canard_radial_cp_offset_m = fixed_spec(0.005);
    opts.performance.engine_thrust_factor = fixed_spec(1);
    opts.performance.act_backlash_deg = fixed_spec(0.25);
    opts.roll_reversal.probability = 1;
    opts.roll_reversal.time_after_takeoff_s = [8, 8];
    opts.roll_reversal.factor_range = [-0.5, -0.5];
end

function spec = fixed_spec(value)
    spec = struct("mean", value, "range", [value, value]);
end

function baseline = geometry_baseline()
    baseline = struct( ...
        "nosecone_pos_x_cp", 0.1, ...
        "body_pos_x_cp", -2, ...
        "tail_pos_x_cp", -5, ...
        "fin_pos_x_cp", -4, ...
        "fin_pos_x_cp_mach2", -4.1, ...
        "canard_chord_root", 0.3, ...
        "canard_height", 0.2, ...
        "canard_sweep_angle", 0, ...
        "canard_spine_width", 0.02, ...
        "canard_spine_effectiveness", 0.5, ...
        "canard_pos_x_roottip", -1, ...
        "rocket_diameter", 0.2, ...
        "canard_radial_offset", 0.02, ...
        "time_idle", 30);
end
