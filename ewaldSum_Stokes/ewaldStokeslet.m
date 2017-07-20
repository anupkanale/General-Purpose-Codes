%% Program to simulate flow past a stokeslet using Ewald Summations
clear; close all

N = 2; % number of stokeslets
rVec = [0.25 0.75;0.5 0.5;0.5 0.5];
fVec = [1 -1; 0 0; 0 0];
alpha = 1;

a1 = [1;0;0];
a2 = [0;1;0];
a3 = [0;0;1];
L = 1;
nReal = 0;
nImag = 5;
u = zeros(N,1);

for mm=1:N
    rVecM = rVec(:,mm);
    
    uReal = realSum(N,rVecM,rVec,fVec,alpha,L,a1,a2,a3,nReal);
%     uFourier = fourierSum(N,rVecM,rVec,fVec,alpha,L,a1,a2,a3,nImag);

%     u(mm) = uReal + uFourier;
end