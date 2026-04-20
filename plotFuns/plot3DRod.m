function name = plot3DRod(g, rout, alpha, color)
% g: frames
% T0: Plot frame, g is considered to be expressed in T0

%% Initialization & Check inputs
if nargin < 3
    alpha = 1;
end

if nargin < 4
    color = [1,0,1];
end

res = size(g,3);

%% Compute bishop frames
% for i = 2:res
%     w = LogSO3(g(1:3,1:3,i-1)'*g(1:3,1:3,i));
%     g(1:3,1:3,i) = g(1:3,1:3,i)*LargeSO3([0,0,-w(3)]);
% end

%% Build meshs
m = 20; % points on cross sectional circle

X = zeros(2*res+1, m);
Y = zeros(2*res+1, m);
Z = zeros(2*res+1, m);
X = zeros(res, m);
Y = zeros(res, m);
Z = zeros(res, m);

t = linspace(0,2*pi,m);
outCircle = [rout*cos(t); rout*sin(t); zeros(1,m);ones(1,m)];
% inCircle = [rin*cos(t); rin*sin(t); zeros(1,m);ones(1,m)];
for i = 1:res
    outCross = g(:,:,i)*outCircle;
    
    X(i,:) = outCross(1,:);
    Y(i,:) = outCross(2,:);
    Z(i,:) = outCross(3,:);

end

name = surf(X,Y,Z, 'MeshStyle', 'both', 'LineStyle', 'none', 'FaceAlpha', alpha, 'FaceColor', color);

end

