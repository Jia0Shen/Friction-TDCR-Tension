function [P1,P2,J1,J2] = getP_J12(tdcr, u, R, p, s)

    % centerline point of P1 P2
    idx1 = round((tdcr.N+1)/2);
    p1c = p(:,idx1);
    p2c = p(:,end);

    R1c = R(:,:,idx1);
    R2c = R(:,:,end);

    offset = [tdcr.r, 0, 0]';

    P1 = p1c + R1c * offset;
    P2 = p2c + R2c * offset;

    % get the jacobian at the centerline
    [Jp, Jw] = computeJacobian(u,R,p,s);

    Jp1 = Jp(3*idx1-2:3*idx1, :);
    Jw1 = Jw(3*idx1-2:3*idx1, :);
    Jp2 = Jp(end-2:end, :);
    Jw2 = Jw(end-2:end, :);

    % find derivative of < P1 = p1c + R1c * offset >;

    J1 = Jp1 - R1c*hat(offset)*Jw1;
    J2 = Jp2 - R2c*hat(offset)*Jw2;

end