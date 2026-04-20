function g = NonLinearLCPFunc(tdcr, prev_state, T0, x, z)
% x = [l0,l1,theta0,theta1,F1x,F1y,F1z,F2x,F2y,F2z,T1,T2]'
% z = [beta1', betat2', lambda1, lambda2]'

l0 = x(1);
l1 = x(2);
theta0 = x(3);
theta1 = x(4);
F1 = x(5:7);
F2 = x(8:10);
T1 = x(11);
T2 = x(12);

beta1 = z(1:2);
beta2 = z(3:4);
lambda1 = z(5);
lambda2 = z(6);

D = [1 -1];
mu = tdcr.mu;

l0_prev = prev_state.l0;
l1_prev = prev_state.l1;

% v0 = (l0-l0_prev) + (l1-l1_prev);
% v1 = l1-l1_prev;

v0 = (l0_prev-l0) + (l1_prev-l1);
v1 = l1_prev-l1;

g = [lambda1*[1;1]+D'*v0;
     lambda2*[1;1]+D'*v1;
     1/2*T0*(exp(mu(1)*theta0)-exp(-mu(1)*theta0)) - [1,1]*beta1;
     1/2*T1*(exp(mu(2)*theta1)-exp(-mu(2)*theta1)) - [1,1]*beta2];


end