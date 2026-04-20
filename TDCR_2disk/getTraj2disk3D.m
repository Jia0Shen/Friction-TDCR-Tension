function [traj_state, traj_stateFl] = getTraj2disk3D(tdcr, tdcr_fl, T0_list)

% generate a traj given T0_list input.


% tdcr_fl = tdcr;
% tdcr_fl.mu = 0*tdcr.mu;

u_ini = zeros(3,tdcr.N);
[TT_ini, R_ini, pc_ini] = solveShape(tdcr.base, u_ini, tdcr.s);
[J_ini, P_ini, l_ini, L_ini] = get_P_J(tdcr, u_ini,R_ini,pc_ini,tdcr.s);

iniState.J = J_ini;
iniState.P = P_ini;
iniState.L = L_ini;
iniState.l = l_ini;
iniState.u = u_ini;
iniState.pc = pc_ini;
iniState.TT = TT_ini;
iniState.x_sol = [P_ini; zeros(5*tdcr.m*tdcr.n,1)];    % x = [P;F;T;theta]
iniState.z_sol = zeros(3*tdcr.m*tdcr.n, 1);
iniState.F = zeros(3*tdcr.m*tdcr.n, 1);
iniState.T = zeros(tdcr.m*tdcr.n, 1);
iniState.theta = zeros(tdcr.m*tdcr.n, 1);
iniState.v = zeros(size(L_ini));

prevState = iniState;
prevState_fl = iniState;

traj_state = {};
traj_stateFl = {};

n_steps = size(T0_list, 2);

for iter = 1:n_steps

    T0_All_i = T0_list(:,iter);

    [updateState] = stepUpdate_t(tdcr, prevState, T0_All_i);
    [updateState_fl] = stepUpdate_t(tdcr_fl, prevState_fl, T0_All_i);
    % updateState_fl = prevState_fl;

    % update
    prevState = updateState;
    prevState_fl = updateState_fl;

    traj_state{end+1} = updateState;
    traj_stateFl{end+1} = updateState_fl;
    
end


end