function [plots] = plot_state(dataset, varargin)
    % plot simulation data on a dashboard for estimator visualization
    
     if nargin == 1 || nargin == 2 || nargin == 3
        tiledlayout(3,2,'TileSpacing','Compact');
        plots.q = nexttile; plots.w = nexttile; plots.v = nexttile; plots.alt = nexttile; 
        plots.cl = nexttile; plots.delta = nexttile;
     elseif nargin == 4
        plots = varargin{3};
     end

    if nargin >= 2
        name = varargin{1};
    else
        name = "";
    end

    time = seconds(dataset.Time);
    
    % names = append(["qw","qx","qy","qz"],name);
    names = ["$q_w$","$q_x$","$q_y$","$q_z$"];
    for i = 1:4
        plot(plots.q, time, dataset.q(:,i), 'DisplayName', names(i));
        hold(plots.q, 'on')
    end
    legend(plots.q, 'show', 'Orientation','vertical', Location='eastoutside' )
    title(plots.q, "Quaternion",'FontWeight','Normal')
    xlabel(plots.q, "Time [s]")
    grid(plots.q, 'on')

    % names = append(["wx","wy","wz"],name);
    names = ["$\omega_x$","$\omega_y$","$\omega_z$"];
    for i = 1:3
        plot(plots.w, time, dataset.w(:,i), 'DisplayName', names(i))
        hold(plots.w, 'on')
    end
    legend(plots.w, 'show', 'Orientation','vertical', Location='eastoutside' )
    title(plots.w, "Angular velocity [rad/s]",'FontWeight','Normal')
    xlabel(plots.w, "Time [s]")
    grid(plots.w, 'on')

    % names = append(["vx","vy","vz"],name);
    names = ["$v_x$","$v_y$","$v_z$"];
    for i = 1:3
        plot(plots.v, time, dataset.v(:,i), 'DisplayName', names(i))
        hold(plots.v, 'on')
    end
    legend(plots.v, 'show', 'Orientation','vertical', Location='eastoutside' )
    title(plots.v, "Velocity [m/s]",'FontWeight','Normal')
    xlabel(plots.v, "Time [s]")
    grid(plots.v, 'on')

    % names = append("alt",name);
    names = "$r_x$";
    for i = 1
        plot(plots.alt, time, dataset.alt(:,i)/1000, 'DisplayName', names(i))
        hold(plots.alt, 'on')
    end
    legend(plots.alt, 'show', 'Orientation','vertical', Location='eastoutside' )
    title(plots.alt, "Altitude [km]",'FontWeight','Normal')
    % legend(plots.alt, 'show', 'Orientation','vertical', Location='eastoutside' )
    title(plots.alt, "Altitude [km]",'FontWeight','Normal')
    xlabel(plots.alt, "Time [s]")
    grid(plots.alt, 'on')

    % names = append("CL",name);
    names = "$C_{L_\delta}$";
    for i = 1
        plot(plots.cl, time, dataset.cl(:,i), 'DisplayName', names(i))
        hold(plots.cl, 'on')
    end
    legend(plots.cl, 'show', 'Orientation','vertical', Location='eastoutside' )
    title(plots.cl, "Canard Coefficient",'FontWeight','Normal')
    xlabel(plots.cl, "Time [s]")
    grid(plots.cl, 'on')

    % names = append("delta",name);
    names = "$\delta$";
    for i = 1
        plot(plots.delta, time, rad2deg(dataset.delta(:,i)), 'DisplayName', names(i))
        hold(plots.delta, 'on')
    end
    legend(plots.delta, 'show', 'Orientation','vertical', Location='eastoutside' )
    title(plots.delta, "Canard Angle [deg]",'FontWeight','Normal')
    xlabel(plots.delta, "Time [s]")
    grid(plots.delta, 'on')

    if  nargin == 3 || nargin == 4
        enablehold = varargin{2};
    else
        enablehold = 'off';
    end
    hold([plots.q, plots.w, plots.v, plots.alt, plots.cl, plots.delta], enablehold)
end