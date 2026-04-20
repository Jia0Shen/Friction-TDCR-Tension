function handles = creatPlotTDCR_3d(tdcr, currState, color, trans)

if nargin < 3
    color = [];
end

if nargin < 4
    trans = 1;
end

    width = 0.5;
    % disk width

    % set colors
    if isempty(color)
        c_disk = [179 163 105]/256;
        c_backbone = [0 48 87]/256;
        c_tendon = [[0 48 87]/256, trans];
        c_circle = [[0 0 0], trans];
    else
        colorrgb = color;
        c_disk = colorrgb;
        c_backbone = colorrgb;
        c_tendon = [colorrgb, trans];
        c_circle = [colorrgb, trans];
    end

    % rr = 1.25*tdcr.r_disk;
    rr = 1.25*tdcr.r;

    TT= currState.TT;
    pc = currState.pc;
    R = TT(1:3,1:3,:);

    x = currState.x_sol;
    z = currState.z_sol;

    [P,F,T,theta,beta,lamb] = stateExtract_t(tdcr, x, z);

    % h1 = plot3Mat(pc, color, 6*width);
    h1 = plot3DRod(TT, 4, trans, c_backbone);

    m = tdcr.m;
    n = tdcr.n;

    h_tendon = [];
    h_disk = [];
    h_hole = [];

    for j = 1:m
        Pj = P(3*n*(j-1)+1:3*n*j);
        Pj_mat = reshape([tdcr.P0(:,j)-[0;0;40];tdcr.P0(:,j); Pj], 3, []);   % add P0
        ht_I = plot3Mat(Pj_mat, c_tendon, 3*width);
        if trans == 1
            ht_I.Marker = '.';
        end
        h_tendon = [h_tendon, ht_I];
    end

    alpha = linspace(0, 2*pi, 100);
    base_circle = [rr*cos(alpha); rr*sin(alpha); zeros(size(alpha))];
    base_holes = tdcr.P0;

    idx1 = ((0:n)*tdcr.NumSD)+1;

    pn = pc(:,idx1);
    Rn = R(:,:,idx1);

    for i = 1:n+1
        pci = pn(:, i);
        Ri = Rn(:,:,i);

        circleI = pci + Ri*base_circle;

        % consider the twist
        if i == 1
            holeI = pci + Ri*base_holes;
        else
            R_twist = RotZ(tdcr.twist_angle(i-1));
            holeI = pci + Ri*R_twist*base_holes;
        end

        hd_I = plot3Mat(circleI, c_circle, 3*width);
        hdp_I = patch(circleI(1,:),circleI(2,:),circleI(3,:), c_disk, 'FaceAlpha', trans);
        % hdp_I = patch(circleI(1,:),circleI(2,:),circleI(3,:), [0 0.4470 0.7410], 'FaceAlpha', 1);
        h_disk = [h_disk, hd_I, hdp_I];

    end


    handles = [h1, h_tendon, h_disk, h_hole];

end