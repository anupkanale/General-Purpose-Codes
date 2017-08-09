% Program to illustrate the use of Ewald summations to calculate total
% electrostatic energy of a system of charges due to Coloumbic interaction
clear; %clc;

a1 = [1;0;0];  % lattice vectors
a2 = [0;1;0];
a3 = [0;0;1];
eps0 = 8.85e-12;
kB = 1.38e-23;

angst = 1e-10;
charge = 0.4238*1.6e-19;
kB4pieps0_Inv = 1/(kB*4*pi*eps0);

layers = 10;
nCharges = 2*(2*layers+1)^3; % number of charged particles
L = 1;

% periodicity vector
pMat = getP(a1,a2,a3,layers);

rVec = zeros(3,nCharges);
q = zeros(nCharges,1);

for ii=1:(2*layers+1)^3
    pVec = pMat(:,ii);
    
    rVec(:,2*ii-1) = [-0.1;0;0] + pVec;
    q(2*ii-1) = -1;
    
    rVec(:,2*ii) = [0.1;0;0] + pVec;
    q(2*ii) = 1;
end

% Visualize charges in space
% for ii=1:nCharges
%     scatter3(rVec(1,ii), rVec(2,ii), rVec(3,ii));
%     hold on;
% end

%% Direct Computation
U = 0;
for ii=1:nCharges
    for jj=ii+1:nCharges
        rij = rVec(:,ii) - rVec(:,jj);
        dr = norm(rij);
        U = U +kB4pieps0_Inv* charge^2/angst* q(ii)*q(jj)/dr;
    end
end

%% Ewald summation
nReal = 10;
nImag = 10;
alpha = 5.6/L; % Ewald parameter

USelf    = -kB4pieps0_Inv* charge^2/angst* sum(q.^2)*alpha/sqrt(pi);
UReal    =  kB4pieps0_Inv* charge^2/angst* realSum(a1,a2,a3,rVec,q,nCharges,nReal,L,alpha);
UFourier =  kB4pieps0_Inv* charge^2/angst* waveSum(a1,a2,a3,rVec,q,nCharges,nImag,L,alpha);

Utot = UReal + UFourier + USelf;