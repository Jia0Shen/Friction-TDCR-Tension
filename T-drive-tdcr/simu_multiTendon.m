clear; clc

tdcr.r = 50;          % disk radii
tdcr.l = 100;          % distance between disks
tdcr.kb = 20e4;

nu = 0.3;
tdcr.kt = tdcr.kb / (1+nu);

% discretized points for backbone
tdcr.m = 2;    % number of tendons
tdcr.n = 5;    % number of disks
tdcr.NumSD = 20;    % devide the section into NumSD pieces of arc-length
tdcr.N = tdcr.n*tdcr.NumSD+1;  % discretized points
tdcr.s = linspace(0, tdcr.n*tdcr.l, tdcr.N);
tdcr.base = eye(4);
tdcr.twist_angle = deg2rad(zeros(1,tdcr.n));

tdcr.Kmat = getKmat(tdcr.kb,tdcr.kt,tdcr.N);
tdcr.inv_K = getKmat(1/tdcr.kb,1/tdcr.kt,tdcr.N);

tdcr.mu = 0.5 * ones(1,tdcr.m*tdcr.n);       % COF

angleBase = (0:tdcr.m-1)/tdcr.m*2*pi;
tdcr.P0 = [tdcr.r*cos(angleBase); tdcr.r*sin(angleBase); zeros(1, tdcr.m)];

% copy another frictinoless tdcr
tdcr_fl = tdcr;
tdcr_fl.mu = 0.0 * ones(1,tdcr.m*tdcr.n);

% input a T0 list
Tmax = 15;
T_mid = Tmax/2;

% 1. linear tension

% T0_list = geneLineTraj(Tmax, 2, 500);
NN = 200;
% T0_list = [zeros(1,NN), linspace(0, Tmax, NN);
%            linspace(0, Tmax, NN), Tmax*ones(1,NN)];
% T0_list = [linspace(0, Tmax, NN), Tmax*ones(1,NN), linspace(Tmax,0,NN), zeros(1,NN); 
%            zeros(1,NN), linspace(0, Tmax, NN), Tmax*ones(1,NN), linspace(Tmax,0,NN)];
% T0_list = square_spiral(Tmax, 2, 0.2)';
dT = 2e-2;
T0_list = LinearInterplNdim([0, 1/2*Tmax, 1/2*Tmax; ...
                             0, 0, Tmax], dT);
% T0_list = LinearInterplNdim([0, 3/4*Tmax, 3/4*Tmax, 0; ...
%                              0, 0, 1/4*Tmax, 0], dT);
% T0_list = linspace(0, Tmax, 1000);
t = linspace(0, 10, size(T0_list,2));

% 2. circle - hysteresis

% t = linspace(0, 10, 200);
% fh = 1;
% tau = 0.4;
% T0_list = T_mid*exp(-tau*t).*sin(2*pi*t*fh-pi/2)+T_mid;

n_steps = size(T0_list, 2);

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

prevState = iniState;
prevState_fl = iniState;

traj_state = {};
traj_stateFl = {};

frictionlessFlag = true;

tic 
for iter = 1:n_steps

    % T0_i = T0_list(iter);
    % T0_All_i = zeros(1,tdcr.m);
    % T0_All_i(1) = T0_i;

    T0_All_i = T0_list(:, iter)';

    % T0_All_i = T0_list(:,iter);

    % frictionpart
    [updateState] = stepUpdate_t(tdcr, prevState, T0_All_i);
    prevState = updateState;
    traj_state{end+1} = updateState;

    % frictionless part
    if frictionlessFlag
        [updateState_fl] = stepUpdate_t(tdcr_fl, prevState_fl, T0_All_i);
        prevState_fl = updateState_fl;
        traj_stateFl{end+1} = updateState_fl;
    end
    
end

time = toc;

disp('Avg CPU time per step is ' + string(time/n_steps))

%% plot 3D / video

plotFun = @creatPlotTDCR_3d;

vid_flag1 = true;

figure()
hold on

hs = plotFun(tdcr, iniState);
plotConfig3D([640 640])

if frictionlessFlag
    hs_fl = plotFun(tdcr_fl, iniState);
end

xlabel('x(mm)');
ylabel('y(mm)');
zlabel('z(mm)');
axis equal

zlim([-60, 20+tdcr.n*tdcr.l])
ylim([-80 80])
xlim([-300 300])

view(25, 10)

if vid_flag1
    vid1 = VideoWriter('out\vid_3disk2tendon.mp4', 'MPEG-4');
    vid1.FrameRate = min(round(n_steps/10), 110);
    open(vid1);
end

for i = 1:n_steps
    iState = traj_state{i};
    delete(hs);
    hs = plotFun(tdcr, iState, [1 0 0], 0.5);

    title(sprintf('T1 = %.3f N, T2 = %.3f N.', T0_list(1, i), T0_list(2, i)))

    if frictionlessFlag
        iState_fl = traj_stateFl{i};
        delete(hs_fl);
        hs_fl = plotFun(tdcr_fl, iState_fl, [0 0 1], 0.5);
    end

    if vid_flag1
        writeVideo(vid1, getframe(gcf));
    else
        pause(0.001)
    end

end

if vid_flag1
    close(vid1);
end

% grid off
% box off
% axis off

%% 2D plots

angle_list = ones(1,n_steps);
angle_listFl = ones(1,n_steps);
beta_list = ones(tdcr.m, n_steps);
beta_listFL = ones(tdcr.m, n_steps);

for i = 1:n_steps
    iState1 = traj_state{i};
    iState_fl = traj_stateFl{i};

    % xSol_i = iState.x_sol;
    % zSol_i = iState.z_sol;

    % END angle
    R_endI = iState1.TT(1:3,1:3,end);
    Angle_endi = acos(R_endI(1));
    angle_list(i) = Angle_endi;

    R_endIFl = iState_fl.TT(1:3,1:3,end);
    angle_listFl(i) = acos(R_endIFl(1));

    % exposed length
    beta_list(:,i) = tdcr.n*tdcr.l - iState1.L(1,:)';
    beta_listFL(:,i) = tdcr.n*tdcr.l - iState_fl.L(1,:)';

end

% figure()
% hold on
% plot(t,T0_list)
% plot(t,T1_list)
% plot(t,T2_list)
% legend('T0', 'T1', 'T2')
% 
% figure()
% hold on
% plot(t, v0_list)
% plot(t, v1_list)
% legend('v0', 'v1')

figure()
hold on
plot(t, T0_list(1,:))
plot(t, T0_list(2,:))

figure()
hold on
plot(t, angle_list/pi*180, 'r')
plot(t, angle_listFl/pi*180, 'b')
xlabel('t (s)')
ylabel('\theta (Deg)')

figure()
hold on
plot(- diff(T0_list, 1), angle_list/pi*180, 'r')
plot(- diff(T0_list, 1), angle_listFl/pi*180, 'b')
xlabel('T0_1 - T0_2 (N)')
ylabel('\theta (Deg)')
xlim([-42 42])
ylim([0 180])
makeVid2D(gca, 'outs\T-drive-outs\vid_T_Bend', 50)

