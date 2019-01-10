% Program showing streamlines past a Stokeslet
%
% Anup Kanale, Nov 25, 2016
%-------------------------------------------------
clc; clear all;
close all;

eta = 1;
eps = 0.01;

rx = linspace(-1,1,21);
ry = linspace(-1,1,21);
F = [1; 0];
for ii=1:length(rx)
    for jj=1:length(ry)
        r = [rx(ii); ry(jj)];
        if norm(r)==0
            r=[eps; eps]; % this is done to avoid the singularity at (0,0)
            continue;
        end
        rDotF = r(1)*F(1) + r(2)*F(2);
        u(ii,jj) = 1/(8*pi*eta*norm(r)) *(F(1) + rDotF*r(1)/norm(r)^2);
        v(ii,jj) = 1/(8*pi*eta*norm(r)) *(F(2) + rDotF*r(2)/norm(r)^2);
    end
end

%Plot streamlines
scale=1000;
set(figure, 'Position', [400, 100, 500, 400]);
plot(rx(11), ry(11), 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'red');
hold on;
q = quiver(rx,ry,u*scale,v*scale);
q.Color = 'blue';
q.LineWidth = 2;

title('Streamlines for slow past a Stokeslet');
axis equal

domain = -1:0.1:1;
[rx,ry] = meshgrid(domain,domain);
starty = -1:0.4:1;
startx=-1*ones(size(starty));
h = streamline(rx,ry,u,v,startx,starty);
set(h, 'LineWidth', 2, 'Color', 'red');
