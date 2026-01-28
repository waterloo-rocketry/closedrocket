%% Configure and Plot
clear 
folderpath = 'monte-carlo/single_uncontrolled/';
name = 'uncon';

%%% Plots

load(append(folderpath, 'result.mat'), "sdt", "sdt_vars");

set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter', 'latex')

tlim = [0, 40, 10];

f_anim = figure(5);
anim = plot_animation(sdt, tlim);

f_sim = figure(1);
plots = plot_state(sdt.rocket_dt, tlim, [2:3]);

% f_est = figure(2);
% plot_state(sdt.est, [0, 40, 10], [1:2, 6]);
% 
% f_err = figure(3);
% plot_state(sdt.error, [0, 50, 0], [1:3]);

% f_cmd = figure(4);
% stairs(sdt.control.Time, rad2deg([sdt.control.(3)]))
% legend("Reference", "Roll angle", "Command")
% ylabel("Angle [deg]")
% title("Command [deg]",'FontWeight','Normal')
% grid on
% exportgraphics(f_cmd, append(folderpath, 'ctrl_cmd.pdf'))

% f_cov = figure(5);
% plot(seconds(sdt.P_norm.Time), sdt.P_norm.(1)(:,1))

fontsize(f_sim, 12, "points")
% fontsize(f_est, 12, "points")
% fontsize(f_err, 12, "points")

% set(f_sim,'units','centimeters','position',[1,1,25,18])
% set(f_est,'units','centimeters','position',[1,1,25,18])
% set(f_err,'units','centimeters','position',[1,1,25,18])


%% export

% exportgraphics(f_sim, append(folderpath, 'sim_', name, '_sim.pdf'), 'ContentType', 'vector')
% exportgraphics(f_est, append(folderpath, 'sim_', name, '_est.pdf'), 'ContentType', 'vector')
% exportgraphics(f_err, append(folderpath, 'sim_', name, '_err.pdf'), 'ContentType', 'vector')

%% animation plots
f_sim;
drawnow
ts = tlim(1):(1/30):tlim(2);
vid = VideoWriter(append(folderpath, 'sim_', name, '_sim.mp4'), 'MPEG-4');
vid.Quality = 95;
vid.FrameRate = 30;
open(vid)
t = 1;
mov_line_w = xline(plots.w, ts(1), 'HandleVisibility','off');
mov_line_v = xline(plots.v, ts(1), 'HandleVisibility','off');
timer = tic;
while t <= length(ts)
    % figure(1)    
    % if toc(timer) >= ts(t)
        mov_line_w.Value = ts(t);
        mov_line_v.Value = ts(t);
        drawnow
        writeVideo(vid, getframe(gcf));
        t = t+1;
    % end
end
close(vid)

%% animation video
%%% Record
f_anim;
anim.VideoRecord = 'on';
anim.VideoQuality = 70;
anim.VideoCompression = 'MPEG-4';
anim.VideoFilename = append(folderpath, 'sim_', name, '_animation.mp4');
anim.VideoTStart = tlim(1);
anim.VideoTFinal = tlim(2);

%%% play
anim.TSTart = tlim(1);
anim.TFinal = tlim(2);
anim.play();


anim.wait();
anim.hide();
anim.VideoRecord = 'off';
anim.delete();

%%
% set(groot, 'defaultAxesTickLabelInterpreter','remove')
% set(groot, 'defaultLegendInterpreter','remove')
% set(groot, 'DefaultTextInterpreter', 'remove')