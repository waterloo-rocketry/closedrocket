function h = plot_animation(dataset, varargin)
    % Plays animation and saves video file
    
    time = seconds(dataset.rocket_dt.Time);
    tlim = [time(1), time(end)];
    if nargin >= 2 && ~isempty(varargin{1})
        time_in = varargin{1};
        tlim = time_in(1:2);
        if length(time_in) == 3
            time = time - time_in(3);
        end
    end

    %%% get Euler angles
    euler = zeros(height(dataset.rocket_dt.q),3);
    for t = 1:height(dataset.rocket_dt.q)
        q = dataset.rocket_dt.q(t,:)';
        euler(t,:) = quaternion_to_euler(q)';
    end
    tt.euler = timetable(dataset.rocket_dt.Time, euler, 'VariableNames', "euler");
    
    % eulerhat = zeros(height(sdt.est.q),3);
    % for t = 1:height(sdt.est.q)
    %     qhat = sdt.est.q(t,:)';
    %     eulerhat(t,:) = quaternion_to_euler(qhat)';
    % end
    % tt.eulerhat = timetable(sdt.rocket_dt.Time, eulerhat, 'VariableNames', "euler");

    
    %%% Animation
    
    h = Aero.Animation;
    h.FramesPerSecond = 30;
    h.TimeScaling = 1;
    h.createBody('testrocket.ac', 'Ac3d');
    % h.createBody('testrocket.ac', 'Ac3d');
    h.createBody('ac3d_xyzisrgb.ac', 'Ac3d');
    animationdata = [time, dataset.rocket_dt.alt, dataset.rocket_dt.pos_yz, euler];
    animationdata = fillmissing(animationdata, 'linear');
    h.Bodies{1}.TimeSeriesSource = animationdata;
    % h.Bodies{2}.TimeSeriesSource = animationdata_hat;
    h.Bodies{2}.TimeSeriesSource = [[0;animationdata(end,1)], repmat([animationdata(1,2:4)+[-3,0,0], deg2rad([0, 180, 0])], 2, 1)];
    h.updateBodies(0);
    % h.Camera.PositionFcn = @staticCameraZoom;
    h.Camera.PositionFcn = @doFirstOrderChaseCameraDynamics;
    h.Camera.AimPoint = h.Bodies{1}.Position;
    h.Camera.UpVector = [1, 0, 0];
    h.Camera.Offset = 4*[-20, 20, 30];
    h.updateCamera(0);
    h.show();

    % %%% Record
    % h.VideoRecord = 'on';
    % h.VideoQuality = 70;
    % h.VideoCompression = 'MPEG-4';
    % h.VideoFilename = 'monte-carlo/animation';
    % h.VideoTStart = tlim(1);
    % h.VideoTFinal = tlim(2);

    %%% play
    % h.TSTart = tlim(1);
    % h.TFinal = tlim(2);
    % h.play();
    
    
    % h.wait();
    % h.hide();
    % h.VideoRecord = 'off';
    % h.delete();
end

function staticCameraZoom(~, Bodies, h)

    if ~isempty(Bodies) && isa(Bodies, 'cell') && isa(Bodies{1},'Aero.Body')
    target     = Bodies{1};
    targetPos  = target.Position;
    else
        targetPos = h.AimPoint; % don't change anything
    end
    
    viewExtent = h.ViewExtent;
    
    %--- Extent of view to render
    
    h.xlim = targetPos(1) + viewExtent;
    h.ylim = targetPos(2) + viewExtent;
    h.zlim = targetPos(3) + viewExtent;
    
    %--- Camera aim point for [x,y,z] aim point
    
    h.AimPoint = targetPos;

    h.Position = h.Offset;

    h.ViewAngle = asin(500/norm(targetPos-h.Offset));
end