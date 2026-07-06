classdef WallStepLogger < matlab.System
    % Logs wall-clock dt between executions of this block.
    % Intended for Simulink diagnostic use only.

    properties (Nontunable)
        Ts = 0.0025                 % expected real-time step, seconds
        MaxSamples = 50000          % enough for 60 s at 400 Hz = 24000
        FileName = "wall_step_log.mat"
    end

    properties (Access = private)
        TimerStart
        K
        SimTimeLog
        WallDtLog
        RatioLog
    end

    methods (Access = protected)

        function setupImpl(obj)
            obj.K = 0;
            obj.TimerStart = tic;

            obj.SimTimeLog = nan(obj.MaxSamples, 1);
            obj.WallDtLog  = nan(obj.MaxSamples, 1);
            obj.RatioLog   = nan(obj.MaxSamples, 1);
        end

        function [wall_dt, ratio] = stepImpl(obj, sim_time)
            wall_dt = toc(obj.TimerStart);
            obj.TimerStart = tic;

            ratio = wall_dt / obj.Ts;

            obj.K = obj.K + 1;

            if obj.K <= obj.MaxSamples
                obj.SimTimeLog(obj.K) = sim_time;
                obj.WallDtLog(obj.K)  = wall_dt;
                obj.RatioLog(obj.K)   = ratio;
            end
        end

        function releaseImpl(obj)
            n = min(obj.K, obj.MaxSamples);

            sim_time = obj.SimTimeLog(1:n);
            wall_dt  = obj.WallDtLog(1:n);
            ratio    = obj.RatioLog(1:n);
            Ts       = obj.Ts;

            save(obj.FileName, "sim_time", "wall_dt", "ratio", "Ts");
        end

        function sts = getSampleTimeImpl(obj)
            sts = createSampleTime(obj, ...
                "Type", "Discrete", ...
                "SampleTime", obj.Ts);
        end

        function [s1, s2] = getOutputSizeImpl(~)
            s1 = [1 1];
            s2 = [1 1];
        end

        function [t1, t2] = getOutputDataTypeImpl(~)
            t1 = "double";
            t2 = "double";
        end

        function [c1, c2] = isOutputComplexImpl(~)
            c1 = false;
            c2 = false;
        end

        function [f1, f2] = isOutputFixedSizeImpl(~)
            f1 = true;
            f2 = true;
        end
    end
end