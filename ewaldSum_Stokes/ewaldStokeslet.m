%% Program to simulate flow past a stokeslet using Ewald Summations
clear; close all

N = 2; % number of stokeslets
nCols = 11;
nRows = 11;
nPoints = nCols*nRows;

lim = 1;
rx = linspace(-lim,lim,nCols);
ry = linspace(-lim,lim,nRows);
rz = linspace(-lim,lim,nZ);
rVec = zeros(3,nPoints);
for ii=1:nRows
    for jj=1:nCols
        pointNum = (ii-1)*nCols+jj;
        rVec(:,pointNum) = [rx(jj); ry(ii); rz(kk)];
    end
end

fVec = zeros(2,nPoints);
% fVec(:,(nPoints+1)/2) = [1; 0];
fLoc = [(nPoints+1)/2 + 2, (nPoints+1)/2 - 2];
fVec(:,fLoc(1)) = [1; 0; 0];
fVec(:,fLoc(2)) = [-1; 0; 0];

a1 = [1;0;0];
a2 = [0;1;0];
a3 = [0;0;1];
L = 1;
nReal = 0;
nImag = 5;
alpha = 0.5;

u = zeros(3,N);
for mm=1:N
    rVecM = rVec(:,mm);

    uReal = realSum(N,rVecM,rVec,fVec,alpha,L,a1,a2,a3,nReal);
    uFourier = fourierSum(N,rVecM,rVec,fVec,alpha,L,a1,a2,a3,nImag);

    u(:,mm) = uReal + uFourier;
end