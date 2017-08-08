%% Program to simulate flow past a stokeslet using Ewald Summations
clear; %close all
tic

% grid points
res = 7;
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
fVec = zeros(3,nPoints);
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
xi = 10;

%% compute flowfield
velVec = zeros(3,nPoints);
uReal  = zeros(3,nPoints);
uFourier = zeros(3,nPoints);
uSelf = zeros(3,nPoints);
for mm=1:nPoints
    rVecM = rVec(:,mm);
    
    uReal(:,mm) = realSum(nStokes,rVecM,r_fVec,fVec,xi,L,a1,a2,a3,nReal);
    uFourier(:,mm) = fourierSum(nStokes,rVecM,r_fVec,fVec,xi,L,a1,a2,a3,nImag);
    uSelf(:,mm) = 4*xi/sqrt(pi) * fVec(:,mm);
    
    velVec(:,mm) = uReal(:,mm) + uFourier(:,mm) - uSelf(:,mm);
end


%% Post-processing
save('ewald_Stokes_data');

% clf(figure(2))
set(figure(), 'position', [1050 1050 850 700]);
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

dummyvel = real(velVec');
dummyr = rVec';
dummyu = u3d';
dummyv = v3d';
dummyreal = uReal';
dummyfourier = uFourier';

toc