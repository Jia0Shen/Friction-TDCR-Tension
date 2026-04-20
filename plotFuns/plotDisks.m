function h_disks = plotDisks(p0_list, R_list, rr, color, trans)

% if nargin

% rr = 1.25*tdcr.r;

alpha = linspace(0, 2*pi, 100);
base_circle = [rr*cos(alpha); rr*sin(alpha); zeros(size(alpha))];

n = size(p0_list, 2);

h_disks = [];

for i = 1:n
    p0i = p0_list(:,i);
    Ri = R_list(:,:,i);

    circleI = p0i + Ri*base_circle;

    hd_I = plot3Mat(circleI, color, 1);
    hdp_I = patch(circleI(1,:),circleI(2,:),circleI(3,:), color, 'FaceAlpha', trans);
    h_disks = [h_disks, hd_I, hdp_I];
end

end