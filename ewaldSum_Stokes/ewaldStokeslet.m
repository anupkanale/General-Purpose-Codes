%% Program to simulate flow past a stokeslet using Ewald Summations
clear; %close all

eta = 1;

% grid points
res = 13;
nX = res;
nY = res*2;
nZ = res;
nPoints = nX*nY*nZ;
L = [1; 1; 1000];
rx = linspace(0,L(1),nX);
ry = linspace(0,L(2),nY);
rz = linspace(0,L(3),nZ);
rVec = zeros(3,nPoints);
for kk=1:nX
    for ii=1:nY
        for jj=1:nZ
            pointNum = (ii-1)*nX + jj + (kk-1)*nX*nY;
            rVec(:,pointNum) = [rx(jj); ry(ii); rz(kk)];
        end
    end
end

% point force location
nStokes = 2; % number of stokeslets
fVec = zeros(3,nPoints);
offset = 5;
fLoc = [(nPoints+1)/2 + offset, (nPoints+1)/2 - offset];
fVec(:,fLoc(1)) = [0.1; 0; 0];
fVec(:,fLoc(2)) = [-0.1; 0; 0];

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
    
    uReal = realSum(nStokes,fLoc,rVecM,rVec,fVec,alpha,L,a1,a2,a3,nReal);
    uFourier = 0;
%     uFourier = fourierSum(nStokes,fLoc,rVecM,rVec,fVec,alpha,L,a1,a2,a3,nImag);
    velVec(:,mm) = uReal + uFourier;
end

%% Post-processing
save('ewald_Stokes_data');
plotGen("ewald");