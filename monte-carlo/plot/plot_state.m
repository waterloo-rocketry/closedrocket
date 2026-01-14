function [plots] = plot_state(dataset, varargin)
    % inputs: dataset, [time limits, time start], plot indices. Plot simulation data on a dashboard for estimator visualization. 
    
    colors(1,:) = [0, 0, 0];  % Black
    colors(2,:) = [0.9, 0.3, 0.1];  % Red
    colors(3,:) = [0.2, 0.8, 0.1];  % Green
    colors(4,:) = [0.1, 0.2, 0.8];  % Blue

    ticknumber = 9;

    time = seconds(dataset.Time);
    tlim = [time(1), time(end)];

    if nargin >= 2 && ~isempty(varargin{1})
        time_in = varargin{1};
        tlim = time_in(1:2);
        if length(time_in) == 3
            time = time - time_in(3);
        end
    end

    if nargin >= 3
        idx = varargin{2};
    else
        idx = 1:6;
    end

    if isscalar(idx)
        tiledlayout(1,1,'TileSpacing','Compact');
    elseif length(idx) == 2
        tiledlayout(2,1,'TileSpacing','Compact');
    elseif length(idx) == 3
        tiledlayout(3,1,'TileSpacing','Compact');
    elseif length(idx) == 4
        tiledlayout(2,2,'TileSpacing','Compact');
    else
        tiledlayout(3,2,'TileSpacing','Compact');
    end 
    
    for i = idx
        switch i
            case 1
                plots.q = nexttile;
            case 2
                plots.w = nexttile;
            case 3
                plots.v = nexttile;
            case 4
                plots.alt = nexttile; 
            case 5
                plots.cl = nexttile;
            case 6
                plots.delta = nexttile;             
        end
    end

    % plots.q = nexttile; plots.w = nexttile; plots.v = nexttile; plots.alt = nexttile; 
    % plots.cl = nexttile; plots.delta = nexttile;
    
    if ismember(1, idx)
    names = ["$q_w$","$q_x$","$q_y$","$q_z$"];
    for i = 1:4
        plot(plots.q, time, dataset.q(:,i), 'DisplayName', names(i), 'Color',colors(i,:));
        hold(plots.q, 'on')
    end
    legend(plots.q, 'show', 'Orientation','vertical', Location='eastoutside', IconColumnWidth=12)
    title(plots.q, "Quaternion",'FontWeight','Normal')
    xlabel(plots.q, "Time [s]")
    xlim(plots.q, tlim)
    % ylim(plots.q,[-1,1])
    xticks(plots.q, linspace(tlim(1), tlim(end), ticknumber) );
    grid(plots.q, 'on')
    plots.q.YAxis.Exponent = 0;
    end


    if ismember(2, idx)
    names = ["$\omega_x$","$\omega_y$","$\omega_z$"];
    for i = 1:3
        plot(plots.w, time, dataset.w(:,i), 'DisplayName', names(i), 'Color',colors(i+1,:))
        hold(plots.w, 'on')
    end
    legend(plots.w, 'show', 'Orientation','vertical', Location='eastoutside', IconColumnWidth=12 )
    title(plots.w, "Angular velocity [rad/s]",'FontWeight','Normal')
    xlabel(plots.w, "Time [s]")
    xlim(plots.w, tlim)
    xticks(plots.w, linspace(tlim(1), tlim(end), ticknumber) );
    grid(plots.w, 'on')
    plots.w.YAxis.Exponent = 0;
    end

    if ismember(3, idx)
    names = ["$v_x$","$v_y$","$v_z$"];
    for i = 1:3
        plot(plots.v, time, dataset.v(:,i), 'DisplayName', names(i), 'Color',colors(i+1,:))
        hold(plots.v, 'on')
    end
    legend(plots.v, 'show', 'Orientation','vertical', Location='eastoutside', IconColumnWidth=12 )
    title(plots.v, "Velocity [m/s]",'FontWeight','Normal')
    xlabel(plots.v, "Time [s]")
    xlim(plots.v, tlim)
    xticks(plots.v, linspace(tlim(1), tlim(end), ticknumber) );
    grid(plots.v, 'on')
    plots.v.YAxis.Exponent = 0;
    % ylim([-200, 800])
    end

    if ismember(4, idx)
    names = "$r_x$";
    for i = 1
        plot(plots.alt, time, dataset.alt(:,i)/1000, 'DisplayName', names(i), 'Color',colors(i+3,:))
        hold(plots.alt, 'on')
    end
    legend(plots.alt, 'show', 'Orientation','vertical', Location='eastoutside', IconColumnWidth=12 )
    title(plots.alt, "Altitude [km]",'FontWeight','Normal')
    % legend(plots.alt, 'show', 'Orientation','vertical', Location='eastoutside' )
    title(plots.alt, "Altitude [km]",'FontWeight','Normal')
    xlabel(plots.alt, "Time [s]")
    xlim(plots.alt, tlim)
    xticks(plots.alt, linspace(tlim(1), tlim(end), ticknumber) );
    grid(plots.alt, 'on')
    plots.alt.YAxis.Exponent = 0;
    % ylim([0, 20])
    end

    if ismember(5, idx)
    names = "$C_{L_\delta}$";
    for i = 1
        plot(plots.cl, time, dataset.cl(:,i), 'DisplayName', names(i), 'Color',colors(i+3,:))
        hold(plots.cl, 'on')
    end
    legend(plots.cl, 'show', 'Orientation','vertical', Location='eastoutside', IconColumnWidth=12 )
    title(plots.cl, "Canard Coefficient",'FontWeight','Normal')
    xlabel(plots.cl, "Time [s]")
    xlim(plots.cl, tlim)
    xticks(plots.cl, linspace(tlim(1), tlim(end), ticknumber) );
    grid(plots.cl, 'on')
    plots.cl.YAxis.Exponent = 0;
    end

    if ismember(6, idx)
    names = "$\delta$";
    for i = 1
        plot(plots.delta, time, rad2deg(dataset.delta(:,i)), 'DisplayName', names(i), 'Color',colors(i+3,:))
        hold(plots.delta, 'on')
    end
    legend(plots.delta, 'show', 'Orientation','vertical', Location='eastoutside', IconColumnWidth=12 )
    title(plots.delta, "Canard Angle [deg]",'FontWeight','Normal')
    xlabel(plots.delta, "Time [s]")
    xlim(plots.delta, tlim)
    xticks(plots.delta, linspace(tlim(1), tlim(end), ticknumber) );
    grid(plots.delta, 'on')
    plots.delta.YAxis.Exponent = 0;
    end
end