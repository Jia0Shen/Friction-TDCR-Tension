function [ff, gg] = NonLinearConstrain_t(tdcr, prevState, T0, x, z)

% p0 = [tdcr.r, 0, 0]';
mu = tdcr.mu;

m = tdcr.m; 
n = tdcr.n;

P0 = tdcr.P0;
ez = [0;0;1];

% ff = zeros(8*m*n,1);
% gg = zeeros(3*m*n, 1);

% f_PJU = zeros(3*m*n, 1);
f_FTP = zeros(3*m*n, 1);
f_TBetaTheta = zeros(m*n, 1);
f_cosTheta = zeros(m*n, 1);
g_lambdaDV = zeros(2*m*n,1);
g_TmuThetaBeta = zeros(m*n,1);

D = [1 -1];
% Kmat = getKmat(tdcr.kb,tdcr.kt,tdcr.N);
inv_K = tdcr.inv_K;

[P,F,T,theta,beta,lmd] = stateExtract_t(tdcr, x, z);

x_prev = prevState.x_sol;
z_prev = prevState.z_sol;
[P_prev,F_prev,T_prev,theta_prev,beta_prev] = stateExtract_t(tdcr, x_prev, z_prev);

l = zeros(m*n,1);
l_prev = zeros(m*n,1);

J = prevState.J;
u_prev = reshape(prevState.u,[],1);

ds = tdcr.s(end) - tdcr.s(end-1);

f_PJU = -P + P_prev - J*u_prev + J*inv_K*J'*F/ds;

% def func for slice
sl = @(x, j, i) x(n*(j-1)+i);
sl2 = @(x, j, i) x(2*(n*(j-1)+i)-1: 2*(n*(j-1)+i));
sl3 = @(x, j, i) x(3*(n*(j-1)+i)-2: 3*(n*(j-1)+i));

% different func
% angle_con = @(theta, t1, t2) cos(theta) + dot(t1/norm(t1), t2/norm(t2));
angle_con = @(theta, t1, t2) pi - theta - atan2(norm(cross(t1,t2)),dot(t1,t2));


for j = 1:m

    % 1. for f

    f_FTP(3*(n*(j-1)+1)-2:3*(n*(j-1)+1)) = - sl3(F,j,1) + sl(T,j,2)*normalize(sl3(P,j,2)-sl3(P,j,1), 'norm') ...
            + sl(T,j,1)*normalize(P0(:,j)-sl3(P,j,1), 'norm');
    f_FTP(3*(n*j)-2:3*(n*j)) = - sl3(F,j,n) + sl(T,j,n)*normalize(sl3(P,j,n-1)-sl3(P,j,n), 'norm');

    % f_cosTheta(n*(j-1)+2) = cos(sl(theta,j,2)) + dot(normalize(sl3(P,j,2)-sl3(P,j,1), 'norm'), normalize(P0(:,j)-sl3(P,j,1), 'norm'));
    % f_cosTheta(n*(j-1)+1) = cos(sl(theta,j,1)) + dot(normalize(sl3(P,j,1)-P0(:,j), 'norm'), -ez);

    f_cosTheta(n*(j-1)+2) = angle_con(sl(theta,j,2), sl3(P,j,2)-sl3(P,j,1), P0(:,j)-sl3(P,j,1));
    f_cosTheta(n*(j-1)+1) = angle_con(sl(theta,j,1), sl3(P,j,1)-P0(:,j), -ez);
    for i = 2:n-1

        f_FTP(3*(n*(j-1)+i)-2:3*(n*(j-1)+i)) = - sl3(F,j,i) + sl(T,j,i+1)*normalize(sl3(P,j,i+1)-sl3(P,j,i), 'norm') ...
            + sl(T,j,i)*normalize(sl3(P,j,i-1)-sl3(P,j,i), 'norm');

        % f_cosTheta(n*(j-1)+i+1) = cos(sl(theta,j,i+1)) + dot(normalize(sl3(P,j,i+1)-sl3(P,j,i), 'norm'), normalize(sl3(P,j,i-1)-sl3(P,j,i), 'norm'));

        f_cosTheta(n*(j-1)+i+1) = angle_con(sl(theta,j,i+1), sl3(P,j,i+1)-sl3(P,j,i), sl3(P,j,i-1)-sl3(P,j,i));
    end


    f_TBetaTheta(n*(j-1)+1) = -sl(T,j,1) + D*sl2(beta,j,1) +1/2*T0(j)*(exp(sl(mu, j, 1)*sl(theta,j,1))+exp(-sl(mu, j, 1)*sl(theta,j,1)));

    for i = 2:n

        f_TBetaTheta(n*(j-1)+i) = -sl(T,j,i) + D*sl2(beta,j,i) +1/2*sl(T,j,i-1)*(exp(sl(mu, j, i)*sl(theta,j,i))+exp(-sl(mu, j, i)*sl(theta,j,i)));

    end

    % find l, l_prev;

    l(n*(j-1)+1) = norm(sl3(P,j,1)-P0(:,j));
    l_prev(n*(j-1)+1) = norm(sl3(P_prev,j,1)-P0(:,j));

    for i = 2:n
        l(n*(j-1)+i) = norm(sl3(P,j,i)-sl3(P,j,i-1));
        l_prev(n*(j-1)+i) = norm(sl3(P_prev,j,i)-sl3(P_prev,j,i-1));
    end

    % 2. for g

    for i = 1:n
        v_j_iMinus1 = sum(l_prev(n*(j-1)+i:n*j)) - sum(l(n*(j-1)+i:n*j));
        g_lambdaDV(2*(n*(j-1)+i)-1: 2*(n*(j-1)+i)) = sl(lmd,j,i)*[1;1] + D'*v_j_iMinus1;
    end

    g_TmuThetaBeta(n*(j-1)+1) = 1/2*T0(j)*(exp(sl(mu, j, 1)*sl(theta,j,1))-exp(-sl(mu, j, 1)*sl(theta,j,1))) - [1 1]*sl2(beta,j,1);
    for i = 2:n
        g_TmuThetaBeta(n*(j-1)+i) = 1/2*sl(T,j,i-1)*(exp(sl(mu, j, i)*sl(theta,j,i))-exp(-sl(mu, j, i)*sl(theta,j,i))) - [1 1]*sl2(beta,j,i);
    end
end


ff = [f_PJU; f_FTP; f_TBetaTheta; f_cosTheta];
gg = [g_lambdaDV; g_TmuThetaBeta];


end

