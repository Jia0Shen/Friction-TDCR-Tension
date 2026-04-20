function [updateState] = solveFrictionStep(tdcr, prevState, T_curr)

    % Kmat = getKmat(tdcr.kb,tdcr.kb,tdcr.n);
    inv_K = tdcr.inv_K;   %diag(1 ./ diag(tdcr.Kmat));
    p0 = [tdcr.r, 0, 0]';

    f = @(x, z) NonLinearConstrain(tdcr, prevState, T_curr, x, z);
    g = @(x, z) NonLinearLCPFunc(tdcr, prevState, T_curr, x, z);

    x_prev = prevState.x_sol;
    z_prev = prevState.z_sol;

    l0_prev = x_prev(1);
    l1_prev = x_prev(2);

    % calculate the dirivative of f,g w.r.t x, z
    [f_x, f_z] = gradientNum2Var(f, x_prev, z_prev);
    [g_x, g_z] = gradientNum2Var(g, x_prev, z_prev);

    % check if f(x_prev, z_prev) = 0

    % f(x_prev, z_prev)

    M = g_z - g_x*inv(f_x)*f_z;
    q = g(x_prev, z_prev) - M * z_prev - g_x*inv(f_x)*f(x_prev,z_prev);

    % solver LCP(M,q).
    [MzPlus_q,z_sol] = LCPSolve(M,q);

    % solve x given z: f(x_sol, z_sol) = 0.

    % x_sol = x_prev + inv(f_x)*f_z*(z_sol-z_prev);   % linearized form (X)
    x_sol = x_prev + inv(f_x)*(-f_z*(z_sol-z_prev)-f(x_prev,z_prev));   % linearized form
    % x_sol = fsolve(@(xx) f(xx, z_sol), x_prev);     % nonlinear form

    % norm(f(x_sol, z_sol))

    % find other states for update next step
    ds = tdcr.s(end) - tdcr.s(end-1);
    
    l0 = x_sol(1);
    l1 = x_sol(2);
    F1F2_sol = x_sol(5:10);
    J_prev = [prevState.J1; prevState.J2];
    u_solList = inv_K*(J_prev'*F1F2_sol)/ds;

    u_sol = reshape(u_solList, 3, []);

    [TT_sol, R_sol, p_sol] = solveShape(tdcr.base, reshape(u_sol,3,tdcr.N), tdcr.s);

    [P1_sol, P2_sol, J1_sol, J2_sol] = getP_J12(tdcr, u_sol, R_sol, p_sol, tdcr.s);

    l1_sol = norm(P2_sol-P1_sol);
    l0_sol = norm(P1_sol - p0);

    theta0_sol = pi - acos((P1_sol-p0)'/norm(P1_sol-p0)*[0;0;-1]);
    theta1_sol = pi - acos((P2_sol-P1_sol)'/norm(P2_sol-P1_sol)*(p0-P1_sol)/norm(p0-P1_sol));

    % update 
    updateState.J1 = J1_sol;
    updateState.J2 = J2_sol;
    updateState.P1 = P1_sol;
    updateState.P2 = P2_sol;
    updateState.l0 = l0;
    updateState.l1 = l1;

    % checkye
    % l0 - l0_sol
    % l1 - l1_sol
    % x_sol(3) - theta0_sol
    % x_sol(4) - theta1_sol
    
    updateState.x_sol = x_sol;
    updateState.z_sol = z_sol;
    updateState.u = u_sol;

    updateState.p = p_sol;
    updateState.TT = TT_sol;

    updateState.v0 = (l0_prev-l0) + (l1_prev-l1);
    updateState.v1 = l1_prev - l1;

    % check the jacobian
    dp1 = updateState.P1 - prevState.P1;
    dp2 = updateState.P2 - prevState.P2;

    F1F2_prev = prevState.x_sol(5:10);
    u_prevList = inv_K*(J_prev'*F1F2_prev);

    u_prev = reshape(prevState.u, [], 1);

    du = u_solList - u_prev;

    dp1_ju = prevState.J1 * du;
    dp2_ju = prevState.J2 * du;

    % err = norm(tdcr.Kmat*u_solList*ds - [J1_sol; J2_sol]'*F1F2_sol)
    % err0 = tdcr.Kmat*u_solList*ds - [J1_sol; J2_sol]*F1F2_sol

    % dp1
    % dp1 - dp1_ju
    % dp2 - dp2_ju

end