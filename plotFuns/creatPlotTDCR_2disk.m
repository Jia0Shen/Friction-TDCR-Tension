function handles = creatPlotTDCR_2disk(tdcr, currState, color, width)

if nargin == 2
    color = 'b';
    width = 0.5;
elseif nargin == 3
    width = 0.5;
end

    % disk width
    rr = 1.25*tdcr.r;

    TT= currState.TT;
    p = currState.p;
    R = TT(1:3,1:3,:);

    x = currState.x_sol;
    l0 = x(1);  l1 = x(2);

    idx1 = round((tdcr.n+1)/2);
    P1c = p(:,idx1);
    P2c = p(:,end);

    R1c = R(:,:,idx1);
    R2c = R(:,:,end);

    offset = [tdcr.r, 0, 0]';
    offset_d = [rr, 0, 0]';

    P0 = offset;
    P1 = P1c + R1c * offset;
    P2 = P2c + R2c * offset;
    
    % points on the disk
    P0d = offset_d;
    P1d = P1c + R1c * offset_d;
    P2d = P2c + R2c * offset_d;

    % the endpoint of the rope
    Pa = [tdcr.r, 0, -(2*tdcr.l-l0-l1)]';

    %
    h1 = plot3Mat(p, color, 6*width);
    h2 = plot3Mat([zeros(3,1), P0d], color, width);
    h3 = plot3Mat([P1c,P1d], color, 4*width);
    h4 = plot3Mat([P2c, P2d], color, 4*width);
    h5 = plot3Mat([P1,P2], color, width);
    h6 = plot3Mat([P0, P1], color, width);
    h7 = plot3Mat([P0, Pa], color, width);

    handles = [h1,h2,h3,h4,h5,h6,h7];

end