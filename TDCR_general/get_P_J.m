function [J, P, l, L] = get_P_J(tdcr, u, R, p, s)

    n = tdcr.n;
    m = tdcr.m;
    kk = tdcr.NumSD;
    twist_angle = tdcr.twist_angle;

    [Jpc3Nx1, Jwc3Nx1] = computeJacobian(u,R,p,s);

    idx1 = ((1:n)*tdcr.NumSD)+1;

    pn = p(:,idx1);
    Rn = R(:,:,idx1);

    P = zeros(3*m*n, 1); 
    J = zeros(3*m*n, 3*tdcr.N); 

    for j = 1:m
        for i = 1:n
            R_twist = RotZ(twist_angle(i));
            % P(3*(n*(j-1)+i)-2:3*(n*(j-1)+i), :) = pn(:, i) + Rn(:,:,i)*tdcr.P0(:,j);
            % J(3*(n*(j-1)+i)-2:3*(n*(j-1)+i), :) = Jpc3Nx1(3*(1+i*kk)-2:3*(1+i*kk), :) - Rn(:,:,i)*hat(tdcr.P0(:,j))*Jwc3Nx1(3*(1+i*kk)-2:3*(1+i*kk), :);
            P(3*(n*(j-1)+i)-2:3*(n*(j-1)+i), :) = pn(:, i) + Rn(:,:,i)*R_twist*tdcr.P0(:,j);
            J(3*(n*(j-1)+i)-2:3*(n*(j-1)+i), :) = Jpc3Nx1(3*(1+i*kk)-2:3*(1+i*kk), :) - Rn(:,:,i)*hat(R_twist*tdcr.P0(:,j))*Jwc3Nx1(3*(1+i*kk)-2:3*(1+i*kk), :);
        end
    end

    % L is the total length vector from P
    l = zeros(n, m);   % l is length of each segments (1 - n)
    L = zeros(n, m);   % L is acumulated length at each disk (0 - n-1)

    for j = 1:m
        Pj = P(1+3*n*(j-1): 3*n*j, :);
        Pj_aug = reshape([tdcr.P0(:,j); Pj], 3, []);
        lj = vecnorm(diff(Pj_aug,1,2), 2, 1);
        l(:, j) = lj';
        L(:, j) = flip(cumsum(flip(lj)))';
    end

end