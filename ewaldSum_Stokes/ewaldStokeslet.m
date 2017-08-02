%% Program to simulate flow past a stokeslet using Ewald Summations
clear; %close all

eta = 1;

% grid points
res = 9;
nX = res;
nY = res;
nZ = 1;
nPoints = nX*nY*nZ;
L = [200000; 200000; 200000];
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
nStokes = 4; % number of stokeslets
fVec = zeros(3,nStokes);
r_fVec = [-0.1 0.1 0 0; 0 0 -0.1 0.1; 0 0 0 0];
fVec(:,1) = [0; 1; 0];
fVec(:,2) = [0; -1; 0];
fVec(:,3) = [1; 0; 0];
fVec(:,4) = [-1; 0; 0];

% ewald summation parameters
a1 = [1;0;0];
a2 = [0;1;0];
a3 = [0;0;1];
nReal = 5;
nImag = 5;
alpha = 0.5;

%% compute flowfield
velVec = zeros(3,nPoints);
for mm=1:nPoints
    rVecM = rVec(:,mm);
    
    uReal = realSum(nStokes,rVecM,r_fVec,fVec,alpha,L,a1,a2,a3,nReal);
    uFourier = 0;
%     uFourier = fourierSum(nStokes,rVecM,r_fVec,fVec,alpha,L,a1,a2,a3,nImag);
    velVec(:,mm) = uReal + uFourier;
end

%% Post-processing
save('ewald_Stokes_data');

clf(figure(2))
set(figure(2), 'position', [1350 0 550 400]);
for ii=1:nStokes
    scatter3(r_fVec(1,ii), r_fVec(2,ii), r_fVec(3,ii), 'filled'); % 'or','markersize', 10, 'markerfacecolor', 'r');
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

% pcolor(rx,ry,velMag2d);
% shading interp;
% colorbar;
scale = 500;
quiv = quiver(rx,ry,u2d*scale,v2d*scale);
quiv.Color = 'blue';
quiv.LineWidth = 1.5;
% axis equal;
grid minor;
xlim([-1 1]); ylim([-1 1]); % zlim([-0.5*L(3), 0.5*L(3)])
view(2);

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