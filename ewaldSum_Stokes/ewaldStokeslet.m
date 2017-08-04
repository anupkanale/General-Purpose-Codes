%% Program to simulate flow past a stokeslet using Ewald Summations
clear; %close all
tic

eta = 1;

% grid points
res = 15;
nX = res;
nY = res*2+1;
nZ = 1;
nPoints = nX*nY*nZ;
L = [2; 2; 2];
rx = linspace(-1, 1, nX);
ry = linspace(-1, 1, nY);
rz = 0;
rVec = zeros(3,nPoints);
for kk=1:nZ
    for jj=1:nY
        for ii=1:nX
            pointNum = ii + (jj-1)*nX + (kk-1)*nX*nY;
            rVec(:,pointNum) = [rx(ii); ry(jj); rz(kk)];
        end
    end
end

% point force location
nStokes = 2; % number of stokeslets
fVec = zeros(3,nStokes);
r_fVec = [-0.1 0.1; 0 0 ; 0 0];
fVec(:,1) = 5*[-1; 0; 0];
fVec(:,2) = 5*[1; 0; 0];
% fVec(:,3) = [1; 0; 0];
% fVec(:,4) = [-1; 0; 0];

% ewald summation parameters
a1 = [1;0;0];
a2 = [0;1;0];
a3 = [0;0;1];
nReal = 10;
nImag = 10;
% alphalist = 5*[1e-3, 1e-2, 1e-1, 1, 1e1, 1e2, 1e3];
% for aa=1:length(alphalist)
alpha = 0.5; %alphalist(ii); %5e-3;

%% compute flowfield
velVec = zeros(3,nPoints);
uReal  = zeros(3,nPoints);
uFourier = zeros(3,nPoints);
parfor mm=1:nPoints
    rVecM = rVec(:,mm);
    
    uReal(:,mm) = realSum(nStokes,rVecM,r_fVec,fVec,alpha,L,a1,a2,a3,nReal);
%     uFourier(:,mm) = 0;
    uFourier(:,mm) = fourierSum(nStokes,rVecM,r_fVec,fVec,alpha,L,a1,a2,a3,nImag);
    velVec(:,mm) = uReal(:,mm) + uFourier(:,mm);
end
dummy4 = uReal';
dummy5 = uFourier';

%% Post-processing
save('ewald_Stokes_data');

% clf(figure(2))
set(figure(2), 'position', [1050 1050 850 700]);
for ii=1:nStokes
    plot(r_fVec(1,ii), r_fVec(2,ii), 'or','markersize', 10, 'markerfacecolor', 'r');
    hold on;
end

% Convert to format appropriate for plotting
u3d = zeros(nX,nY,nZ);
v3d = zeros(nX,nY,nZ);
w3d = zeros(nX,nY,nZ);
velMag3d = zeros(nX,nY,nZ);
for pointNum=1:nPoints
    kk = floor((pointNum-1)/(nX*nY)) + 1;
    jj = floor((pointNum-(kk-1)*nX*nY-1)/nX) + 1;
    ii = rem(pointNum-1,nX) + 1;
    
    u3d(ii,jj,kk) = velVec(1,pointNum);
    v3d(ii,jj,kk) = velVec(2,pointNum);
    w3d(ii,jj,kk) = velVec(3,pointNum);
    velMag3d(ii,jj,kk) = norm(velVec(:,pointNum));
end

% Plot flowfield in 2D
u2d = squeeze(u3d(:,:,(nZ+1)/2));
v2d = squeeze(v3d(:,:,(nZ+1)/2));
velMag2d = squeeze(velMag3d(:,:,(nZ+1)/2));

[x,y] = meshgrid(rx,ry);
% pcolor(x,y,velMag2d');
% shading interp;
% colorbar;
% caxis([0 5]);

quiv = quiver(rx,ry,u3d',v3d');
quiv.Color = 'red';
quiv.LineWidth = 1.5;
% axis equal;
grid minor;
xlim([-1.01 1.01]); ylim([-1.01 1.01]); % zlim([-0.5*L(3), 0.5*L(3)])
% view(2);
title('Ewald method');

dummy = velVec';
dummy2 = u3d';
dummy3 = v3d';

% probe(aa) = v3d(3,5);
% end
% plot(alphalist, probe, 'linewidth', 2)

toc

% 3D plotting
% [x,y,z] = meshgrid(rx,ry,rz);
% quiv = quiver3(x,y,z,u3d,v3d,w3d);
% quiv.Color = 'White';
% quiv.LineWidth = 1.0;
% Slice plot
% xslice = 1;
% yslice = 1;
% zslice = 0.5;
% slice(x,y,z,velMag3d,xslice,yslice,zslice)