function list_xy = square_spiral(Tmax, interval, ds)

% interval is the gap between spiral,
% ds is the step length
    n = Tmax / (2*interval);

    xy = [Tmax, 0];

    list_xy = [0,0; (ds:ds:Tmax)', 0*(ds:ds:Tmax)'];

    for i = 1:n

        xy_prev = xy;
        xy = xy + [0, 40-(2*i-2)*interval];
        dx = ds / (40-(2*i-2)*interval);
        scale = (dx:dx:1)';
        add_list = xy_prev + [scale scale] .* (xy - xy_prev);
        list_xy = [list_xy; add_list];

        xy_prev = xy;
        xy = xy + [-(40-(2*i-2)*interval), 0];
        dx = ds / (40-(2*i-2)*interval);
        scale = (dx:dx:1)';
        add_list = xy_prev + [scale scale] .* (xy - xy_prev);
        list_xy = [list_xy; add_list];
        
        xy_prev = xy;
        xy = xy + [0, -(40-(2*i-1)*interval)];
        dx = ds / (40-(2*i-1)*interval);
        scale = (dx:dx:1)';
        add_list = xy_prev + [scale scale] .* (xy - xy_prev);
        list_xy = [list_xy; add_list];
        
        xy_prev = xy;
        xy = xy + [40-(2*i-1)*interval, 0];
        dx = ds / (40-(2*i-1)*interval);
        scale = (dx:dx:1)';
        add_list = xy_prev + [scale scale] .* (xy - xy_prev);
        list_xy = [list_xy; add_list];


    end
end