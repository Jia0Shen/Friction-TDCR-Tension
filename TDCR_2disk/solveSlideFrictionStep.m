function [updateState] = solveSlideFrictionStep(tdcr, prevState, T_curr)
% This code implements a comparison to the LCP-based method: 
% it only considers the sliding friction.

    % Kmat = getKmat(tdcr.kb,tdcr.kb,tdcr.n);
    inv_K = tdcr.inv_K;   %diag(1 ./ diag(tdcr.Kmat));
    p0 = [tdcr.r, 0, 0]';

    f = @(x) NonLinearConstrainSlide(tdcr, prevState, T_curr, x);

    x_prev = prevState.x_sol;

    l0_prev = x_prev(1);
    l1_prev = x_prev(2);

    % only need solve f
    opt = optimoptions("fsolve", 'Display', 'none');
    x_sol = fsolve(f, x_prev, opt);     % nonlinear form

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
    updateState.T0 = T_curr;
    updateState.J1 = J1_sol;
    updateState.J2 = J2_sol;
    updateState.P1 = P1_sol;
    updateState.P2 = P2_sol;
    updateState.l0 = l0;
    updateState.l1 = l1;
    
    updateState.x_sol = x_sol;
    updateState.u = u_sol;

    updateState.p = p_sol;
    updateState.TT = TT_sol;

    updateState.v0 = (l0_prev-l0) + (l1_prev-l1);
    updateState.v1 = l1_prev - l1;

    % check the jacobian
    % dp1 = updateState.P1 - prevState.P1;
    % dp2 = updateState.P2 - prevState.P2;
    % 
    % F1F2_prev = prevState.x_sol(5:10);
    % u_prevList = inv_K*(J_prev'*F1F2_prev);
    % 
    % u_prev = reshape(prevState.u, [], 1);
    % 
    % du = u_solList - u_prev;
    % 
    % dp1_ju = prevState.J1 * du;
    % dp2_ju = prevState.J2 * du;


end