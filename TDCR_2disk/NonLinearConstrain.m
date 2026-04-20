function y = NonLinearConstrain(tdcr, prev_state, T0, x, z)
% current state (x): [l0,l1,theta0,theta1,F1x,F1y,F1z,F2x,F2y,F2z,T1,T2]'
% z = [beta1', betat2', lambda1, lambda2]'

p0 = [tdcr.r, 0, 0]';
mu = tdcr.mu;

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

ds = tdcr.s(end) - tdcr.s(end-1);

P1P2_prev = [prev_state.P1; prev_state.P2];
u_prev = reshape(prev_state.u, [], 1);
J = [prev_state.J1; prev_state.J2];

D = [1 -1];
% Kmat = getKmat(tdcr.kb,tdcr.kb,tdcr.n);
inv_K = tdcr.inv_K;  % diag(1 ./ diag(tdcr.Kmat));

P1P2 = [p0 + l0*[sin(theta0); 0; cos(theta0)];
     p0 + l0*[sin(theta0); 0; cos(theta0)] + l1*[sin(theta0+theta1); 0; cos(theta0+theta1)]];

pickxz = @(xx) xx([1,3,4,6],:);

y = [pickxz( P1P2 - P1P2_prev - J*inv_K*J'/ds*[F1;F2] + J*u_prev );
     F1 + T1*[sin(theta0); 0; cos(theta0)] - T2*[sin(theta0+theta1); 0; cos(theta0+theta1)];
     F2 + T2*[sin(theta0+theta1); 0; cos(theta0+theta1)];
     T1 - ( 1/2*T0*(exp(mu(1)*theta0)+exp(-mu(1)*theta0))+D*beta1 );
     T2 - ( 1/2*T1*(exp(mu(2)*theta1)+exp(-mu(2)*theta1))+D*beta2 )];


end