clc;
clear;

% load experimental results:
% (force_pull): the input tendon's tension
% (ptop_pull, pmid_pull, RTop_pull, RMid_pull): 
%       the position/orientation of the top/middle disk.
% (BendingExp_pull): bending angle of the top disk
% (timeCam):  time frame

load("out\exp_data.mat")

% if we redo the simulation, 0/1.
retrain = 1;
%% simulation of 3D

E_steel = 53e9;   % 56 GPa
D_tube = 2.3822;
d_tube = 1.9304;  % may not be accurate
I_tube = 1/64*pi*(D_tube^4-d_tube^4);
tdcr.r_disk = 40;
tdcr.l = 406.7/2;
tdcr.kb = E_steel*(1e-6)*I_tube;
nu = 0.3;
tdcr.kt = tdcr.kb / (1+nu);

% discretized points for backbone
tdcr.m = 1;         % 1 disk
tdcr.n = 2;
tdcr.NumSD = 100;    % devide the section into NumSD pieces of arc-length
tdcr.N = tdcr.n*tdcr.NumSD+1;  % discretized points
tdcr.s = linspace(0, tdcr.n*tdcr.l, tdcr.N);
tdcr.base = eye(4);
tdcr.twist_angle = deg2rad([0 0]);

tdcr.K = diag([tdcr.kb, tdcr.kb, tdcr.kt]);
tdcr.Kmat = getKmat(tdcr.kb,tdcr.kt,tdcr.N);
tdcr.inv_K = getKmat(1/tdcr.kb,1/tdcr.kt,tdcr.N);

angleBase = (0:tdcr.m-1)/tdcr.m*2*pi;
tdcr.P0 = [tdcr.r_disk*cos(angleBase); tdcr.r_disk*sin(angleBase); zeros(1, tdcr.m)];

tdcr.nj = tdcr.n;
tdcr.r = {[tdcr.P0, tdcr.P0, tdcr.P0]};
tdcr.oriTendonLen = 2*tdcr.l;

% mus
calib_mu1 = [0.23 0.26];

% frictionless
tdcr_fl = tdcr;
tdcr_fl.mu = [0 0];
tdcr_fl = changeStiffness(tdcr_fl, 59e9*(1e-6)*I_tube);

% discretized tension step length
dT = 500e-3;

% find the T0 list
[T0_pull, Tidx_pull] = LinearInterpl([0;force_pull]*9.8, dT);

% exclude zero
TidxExZero_pull = Tidx_pull(2:end);

if retrain
    
    tdcr.mu = calib_mu1;
    traj_pullSF = getTraj2disk3DSlide(tdcr, T0_pull);
    [traj_pull, traj_pullFl] = getTraj2disk3D(tdcr, tdcr_fl, T0_pull);

    save('out\simu3D_pull.mat')
else
    load('out\simu3D_pull.mat')
end


%% make plots: pull-release 2d
% --- this codes can be stand alone with user-defined "xx_used"
T_used = T0_pull;
traj_used = {traj_pull, traj_pullSF, traj_pullFl};
idx_used = TidxExZero_pull;
exp_used = {ptop_pull, RTop_pull, pmid_pull, RMid_pull, BendingExp_pull};
timeCam_used = timeCam;
force_used = force_pull;
% ---       ---

FrameRate_exp = 5;

traj_lcp = traj_used{1};
traj_sf = traj_used{2};
traj_fl = traj_used{3};

exp_ptop = exp_used{1};
exp_Rtop = exp_used{2};
exp_pmid = exp_used{3};
exp_Rmid = exp_used{4};
exp_Bending = exp_used{5};

n_used = length(T_used);

Bending_pull = zeros(1,n_used);
Bending_pullSf = zeros(1,n_used);
Bending_pullFl = zeros(1,n_used);

tip_pull = zeros(3,n_used);
tip_pullSf = zeros(3,n_used);
tip_pullFl = zeros(3,n_used);

v_pull = zeros(2, n_used);
v_pullSf = zeros(2, n_used);
v_pullFl = zeros(2, n_used);

for i = 1:n_used
    iState = traj_lcp{i};
    iState_sf = traj_sf{i};
    iState_fl = traj_fl{i};

    % END angle
    % tip_pull(:,i) = iState.TT(1:3,4,end);
    [tip_pull(:,i), Bending_pull(i)] = getAngleAndTip(iState.TT(:,:,end));

    [tip_pullFl(:,i), Bending_pullFl(i)] = getAngleAndTip(iState_fl.TT(:,:,end));

    [tip_pullSf(:,i), Bending_pullSf(i)] = getAngleAndTip(iState_sf.TT(:,:,end));

    if i > 1
        v_pull(:, i) = traj_lcp{i-1}.L - traj_lcp{i}.L;
        v_pullFl(:, i) = traj_fl{i-1}.L - traj_fl{i}.L;
        v_pullSf(:, i) = traj_sf{i}.v{1}';
    end
end

% calculate the relative error, our method
bending_avgError = mean(abs(Bending_pull(idx_used) - exp_Bending'));
bending_maxError = max(abs(Bending_pull(idx_used) - exp_Bending'));
tip_avgError = mean(vecnorm(tip_pull(:,idx_used)-exp_ptop, 2, 1));
tip_maxError = max(vecnorm(tip_pull(:,idx_used)-exp_ptop, 2, 1));

% calculate the relative error, sliding friction
bending_avgErrorSf = mean(abs(Bending_pullSf(idx_used) - exp_Bending'));
bending_maxErrorSf = max(abs(Bending_pullSf(idx_used) - exp_Bending'));
tip_avgErrorSf = mean(vecnorm(tip_pullSf(:,idx_used) - exp_ptop, 2, 1));
tip_maxErrorSf = max(vecnorm(tip_pullSf(:,idx_used) - exp_ptop, 2, 1));

% calculate the relative error, Frictionless
bending_avgErrorFl = mean(abs(Bending_pullFl(idx_used) - exp_Bending'));
bending_maxErrorFl = max(abs(Bending_pullFl(idx_used) - exp_Bending'));
tip_avgErrorFl = mean(vecnorm(tip_pullFl(:,idx_used) - exp_ptop, 2, 1));
tip_maxErrorFl = max(vecnorm(tip_pullFl(:,idx_used) - exp_ptop, 2, 1));

disp('---- Pull-release experiment: ------- ')

fprintf('LCP, bending error: Avg. %.4f, Max. %.4f;   Tip error: Avg. %.4f, Max. %.4f \n', ...
    [bending_avgError, bending_maxError, tip_avgError, tip_maxError]);
fprintf('Sliding friction, bending error: Avg. %.4f, Max. %.4f;   Tip error: Avg. %.4f, Max. %.4f \n', ...
    [bending_avgErrorSf, bending_maxErrorSf, tip_avgErrorSf, tip_maxErrorSf]);
fprintf('Frictionless, bending error: Avg. %.4f, Max. %.4f;   Tip error: Avg. %.4f, Max. %.4f \n', ...
    [bending_avgErrorFl, bending_maxErrorFl, tip_avgErrorFl, tip_maxErrorFl]);

disp('------------')

figure()
hold on
fl_select = 1:25;
plot(T_used(idx_used), Bending_pull(idx_used), '-r.', 'LineWidth',1)
plot(T_used(idx_used), Bending_pullSf(idx_used), '-g.', 'LineWidth',1)
plot(T_used(idx_used), Bending_pullFl(idx_used), '-b.', 'LineWidth',1)
plot(force_used*9.8, exp_Bending, '-k.', 'LineWidth',1);
xlabel('T_0 (N)')
ylabel('Bending Angle (deg)')
xlim([0 8.75])
ylim([0 300])

plotConfig2D([400 400])

% if you want to save fig
% saveas(gcf, "out/fig_bend_pull", 'png')
% makeVid2D(gca, "out/vid_bend_pull.mp4", FrameRate_exp)

t_total_pull = max(timeCam_used) - min(timeCam_used);
t_pull = linspace(0, t_total_pull, n_used);
dt_pull = n_used / t_total_pull;

timeCam1From0 = timeCam_used - timeCam_used(1);

figure()
hold on
plot(timeCam1From0, v_pull(1,idx_used)*dt_pull, 'r', 'LineStyle','-', 'LineWidth', 1)
plot(timeCam1From0, v_pull(2,idx_used)*dt_pull, 'r', 'LineStyle',':', 'LineWidth', 2)
plot(timeCam1From0, v_pullSf(1,idx_used)*dt_pull, 'g', 'LineStyle','-', 'LineWidth', 1)
plot(timeCam1From0, v_pullSf(2,idx_used)*dt_pull, 'g', 'LineStyle',':', 'LineWidth', 2)
plot(timeCam1From0, v_pullFl(1,idx_used)*dt_pull, 'b', 'LineStyle','-', 'LineWidth', 1)
plot(timeCam1From0, v_pullFl(2,idx_used)*dt_pull, 'b', 'LineStyle',':', 'LineWidth', 2)
xlim([-0.5, t_total_pull+0.5])
ylim([-36.5, 15])
xlabel('t (s)')
ylabel('v (mm/s)')
plotConfig2D([400 250])

% if you want to save fig
% saveas(gcf, "out/fig_v_pull", 'png')
% makeVid2D(gca, "out/vid_v_pull.mp4", FrameRate_exp)


%% make 3D / video

% vid_flag1 = true;
% 
% figure()
% hold on
% 
% hs = [];
% hs_fl = [];
% hs_sf = [];
% htop = [];
% 
% xlabel('x(mm)');
% ylabel('y(mm)');
% zlabel('z(mm)');
% axis equal
% 
% ylim([-60 60])
% xlim([-50 350])
% zlim([-50, 30+tdcr.n*tdcr.l])
% 
% view(15, 20)
% plotConfig3D([560 560])
% 
% if vid_flag1
%     vid1 = VideoWriter('out\vid3D_pull.mp4', 'MPEG-4');
%     vid1.FrameRate = FrameRate_exp;
%     open(vid1);
% end
% 
% for i = 1:length(idx_used)
%     idx = idx_used(i);
%     iState = traj_lcp{idx};
%     iState_fl = traj_fl{idx};
%     iState_sf = traj_sf{idx};
% 
%     delete(hs);
%     delete(hs_fl);
%     delete(hs_sf);
%     delete(htop);
% 
%     hs = creatPlotTDCR_3d(tdcr, iState, [1 0 0], 0.6);
%     hs_fl = creatPlotTDCR_3d(tdcr, iState_fl, [0 0 1], 0.6);
% 
%     % hs_sf = creatPlotTDCR_3d(tdcr, iState_sf, [0 0 1], 0.6);
% 
%     % plot the orientation
%     Rexp_i = cat(3, exp_Rtop(:,:,i), exp_Rmid(:,:,i), eye(3));
%     pexp_i = [exp_ptop(:, i), exp_pmid(:, i), zeros(3,1)];
% 
%     % hexp = plotCoord(ptop_i, Rtop_i, 40, 'k');
%     hexp = plotDisks(pexp_i, Rexp_i, 1.25*tdcr.r, 'k', 0.7);
% 
%     htop = plotCoord(ptop_i, Rtop_i, 40, 'k');
% 
%     title('T_0 = '+ string(T_used(idx_used(i))) + ' N.')
% 
%     if vid_flag1
%         writeVideo(vid1, getframe(gcf));
%     else
%         pause(0.001)
%     end
% 
% end
% 
% if vid_flag1
%     close(vid1);
% end

%%
function [p_end, angle] = getAngleAndTip(T_end)
    R_end = T_end(1:3, 1:3);
    p_end = T_end(1:3, 4);

    Angle_end = atan2(R_end(1,3), R_end(3,3));
    Angle_end_pi = wrapTo2Pi(Angle_end);
    angle = rad2deg(Angle_end_pi);
end