function [updateState] = stepUpdate_t(tdcr, prevState, T_curr)

    % Kmat = getKmat(tdcr.kb,tdcr.kt,tdcr.N);

    function f_out = f(x, z)
        [f_out, ~] = NonLinearConstrain_t(tdcr, prevState, T_curr, x, z);
    end

    function g_out = g(x, z)
        [~, g_out] = NonLinearConstrain_t(tdcr, prevState, T_curr, x, z);
    end

    x_prev = prevState.x_sol;
    z_prev = prevState.z_sol;

    % calculate the dirivative of f,g w.r.t x, z
    [f_x, f_z] = gradientNum2Var(@f, x_prev, z_prev);
    [g_x, g_z] = gradientNum2Var(@g, x_prev, z_prev);

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

    ds = tdcr.s(end) - tdcr.s(end-1);
    
    J_prev = prevState.J;
    [P_sol,F_sol,T_sol,theta_sol,beta_sol] = stateExtract_t(tdcr, x_sol, z_sol);
    u_solList = tdcr.inv_K*(J_prev'*F_sol)/ds;
    u_sol = reshape(u_solList, 3, []);

    [TT_sol, R_sol, pc_sol] = solveShape(tdcr.base, u_sol, tdcr.s);

    updateState.x_sol = x_sol;
    updateState.z_sol = z_sol;
    updateState.u = u_sol;

    updateState.pc = pc_sol;
    updateState.TT = TT_sol;

    % working on

    [J, P_check, l, L] = get_P_J(tdcr,u_sol,R_sol,pc_sol,tdcr.s);

    % if max(P_check - P_sol) > 1e-3
    %     warning('wrong p_check')
    % end

    updateState.P = P_check;   % no using P_sol
    updateState.F = F_sol;
    updateState.J = J;
    updateState.L = L;
    updateState.l = l;
    updateState.T = T_sol;

    % needs to be updated
    updateState.theta = theta_sol;

    % updateState.v = nan;
    % updateState.l = nan;

end

