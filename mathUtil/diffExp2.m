function [dexpu, C_halfv] = diffExp2(u,ds)

norm_u = norm(u);
% skew_v = so3(v);
% skew_v = [0,-ds,0; ds,0,0; 0,0,0];
skew_v = zeros(3);
skew_v([2,4]) = [ds, -ds];
if norm_u < 1e-12
    dexpu = eye(3);
    C_halfv = 0.5*skew_v;
else
    norm_u_half = 0.5*norm_u;
    norm_u_sqr = (norm_u^2);
    one_over_norm_u_sqr = 1/norm_u_sqr;
    
    s = sin(norm_u_half)/norm_u_half;
    c = cos(norm_u_half);
    
    alpha = s*c;
    beta = s^2;
    
    skew_u = [0,-u(3),u(2); u(3),0,-u(1); -u(2),u(1),0];
    skew_u_sqr = -norm_u_sqr*eye(3) + u*u';
%     skew_u_sqr = [-u(2)^2-u(3)^2, u(1)*u(2), u(1)*u(3);
%                   u(1)*u(2), -u(1)^2-u(3)^2, u(2)*u(3);
%                   u(1)*u(3), u(2)*u(3), -u(1)^2-u(2)^2];
%  
    const1 = one_over_norm_u_sqr*(1-alpha);
    const2 = 3*(1-alpha)*one_over_norm_u_sqr;
    
    dexpu = eye(3) + 0.5*beta*skew_u + const1*skew_u_sqr;
  
    C_halfv = 0.5*(beta)*skew_v + const1 * (skew_v*skew_u + skew_u*skew_v)...
              + (alpha-beta)*one_over_norm_u_sqr*u(3)*ds*skew_u ...
              + one_over_norm_u_sqr*(0.5*beta - const2)*u(3)*ds*skew_u_sqr;
    
end

end

