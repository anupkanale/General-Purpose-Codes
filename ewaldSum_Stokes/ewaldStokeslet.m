%% Program to simulate flow past a stokeslet using Ewald Summations
clear; close all

eta = 1;
% N = 2; % number of stokeslets

% grid points
res = 5;
nX = res;
nY = res;
nZ = res;
nPoints = nX*nY*nZ;
lim = 1;
rx = linspace(-lim,lim,nX);
ry = linspace(-lim,lim,nY);
rz = linspace(-lim,lim,nZ);
rVec = zeros(3,nPoints);
for kk=1:nZ
    for ii=1:nX
        for jj=1:nY
            pointNum = (ii-1)*nX + jj + (kk-1)*nX*nY;
            rVec(:,pointNum) = [rx(jj); ry(ii); rz(kk)];
        end
    end
end

% point force location
fVec = zeros(3,nPoints);
fLoc = [(nPoints+1)/2 + 1, (nPoints+1)/2 - 1];
fVec(:,fLoc(1)) = [1; 0; 0];
fVec(:,fLoc(2)) = [-1; 0; 0];

% Ewald summation problems
a1 = [1;0;0];
a2 = [0;1;0];
a3 = [0;0;1];
L = 1;
nReal = 0;
nImag = 5;
alpha = 0.5;

velVec = zeros(3,nPoints);
for mm=1:nPoints
    rVecM = rVec(:,mm);

    uReal = realSum(nPoints,rVecM,rVec,fVec,alpha,L,a1,a2,a3,nReal);
    uFourier = fourierSum(nPoints,rVecM,rVec,fVec,alpha,L,a1,a2,a3,nImag);

    velVec(:,mm) = uReal + uFourier;
end

%% Plotting
%----------------
scale = 500;
set(figure, 'Position', [2700, 1000, 1000, 700]);
title('Streamlines for slow past a Stokeslet');

for ii=1:length(fLoc)
    scatter3(rVec(1,fLoc(ii)), rVec(2,fLoc(ii)), rVec(3,fLoc(ii)), 'filled');
    hold on;
end
xlim([-1 1]); ylim([-1 1]); zlim([-1 1])

% Convert to format appropriate for plotting
u3d = zeros(nX,nY,nZ);
v3d = zeros(nX,nY,nZ);
w3d = zeros(nX,nY,nZ);
for pointNum=1:nPoints
    kk = floor((pointNum-1)/(nX*nY)) + 1;
    jj = rem(pointNum-1,nX) + 1;
    ii = floor((pointNum-(kk-1)*nX*nY-1)/nX) + 1;
    u3d(ii,jj,kk) = velVec(1,pointNum);
    v3d(ii,jj,kk) = velVec(2,pointNum);
    w3d(ii,jj,kk) = velVec(3,pointNum);
end

% Plot flowfield in 2D
u2d = squeeze(u3d(:,:,(nZ+1)/2));
v2d = squeeze(v3d(:,:,(nZ+1)/2));
quiv = quiver(rx,ry,u2d*scale,v2d*scale);
quiv.Color = 'blue';
quiv.LineWidth = 1.5;
view(2);
axis equal;
grid minor;

save('ewald_Stokes_data')