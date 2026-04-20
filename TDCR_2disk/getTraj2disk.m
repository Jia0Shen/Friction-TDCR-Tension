function [traj_state, traj_state_sf, traj_stateFl] = getTraj2disk(tdcr, tdcr_fl, T0_list)

% generate a traj given T0_list input.
% tdcr.base = g0;

n_steps = size(T0_list, 2);

%initialize states
u_ini = zeros(3,tdcr.N);
[TT_ini, R_ini, p_ini] = solveShape(tdcr.base, u_ini, tdcr.s);
[P1_ini, P2_ini, J1_ini, J2_ini] = getP_J12(tdcr, u_ini, R_ini, p_ini, tdcr.s);

iniState.T0 = 0;

iniState.J1 = J1_ini;
iniState.J2 = J2_ini;
iniState.P1 = P1_ini;
iniState.P2 = P2_ini;
iniState.l0 = tdcr.l;
iniState.l1 = tdcr.l;
iniState.x_sol = [tdcr.l;tdcr.l;zeros(10,1)];
iniState.z_sol = [zeros(4,1); 1; 1];
iniState.u = u_ini;
iniState.p = p_ini;
iniState.TT = TT_ini;
iniState.v0 = 0;
iniState.v1 = 0;

prevState = iniState;
prevState_sl = iniState;
prevState_fl = iniState;

traj_state = {};
traj_state_sf = {};
traj_stateFl = {};

% tdcr_fl = tdcr;
% tdcr_fl.mu = 0;

for iter = 1:n_steps

    T0_i = T0_list(iter);

    [updateState] = solveFrictionStep(tdcr, prevState, T0_i);
    prevState = updateState;
    traj_state{end+1} = updateState;

    % new: add sliding friction
    [updateState_sf] = solveSlideFrictionStep(tdcr, prevState_sl, T0_i);
    prevState_sl = updateState_sf;
    traj_state_sf{end+1} = updateState_sf;

    [updateState_fl] = solveFrictionStep(tdcr_fl, prevState_fl, T0_i);
    prevState_fl = updateState_fl;
    traj_stateFl{end+1} = updateState_fl;

end

end