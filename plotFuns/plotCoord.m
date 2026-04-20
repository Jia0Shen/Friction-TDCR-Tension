function hs = plotCoord(p0, R, ar, color)

% if nargin
    h1 = plot3Mat([p0, p0+ar*R(:,1)], color, 1.5);   % x axis
    h2 = plot3Mat([p0, p0+ar*R(:,2)], color, 1.5);   % y axis
    h3 = plot3Mat([p0, p0+ar*R(:,3)], color, 1.5);   % z axis
    h4 = plot3(p0(1), p0(2), p0(3),'o', 'Color', color);

    hs = [h1 h2 h3 h4];

end