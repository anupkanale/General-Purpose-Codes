% Program to illustrate the use of Ewald summations to calculate total
% electrostatic energy of a system of charges due to Coloumbic interaction
tic
clear; %clc;

a1 = [1;0;0];  % lattice vectors
a2 = [0;1;0];
a3 = [0;0;1];
eps0 = 8.85e-12;
kB = 1.38e-23;

angst = 1e-10;
charge = 0.4238*1.6e-19;
kB4pieps0_Inv = 1/(kB*4*pi*eps0);

layers = 20;
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

%% Direct Computation
set(figure(200), 'position', [50 50 1000 800]);
nList = [1,5,7];
Udirect = zeros(length(nList),1);
parfor kk=1:length(nList)
    Usum = 0;
    for ii=1:nCharges
        ri = rVec(:,ii);
        qi = q(ii);
        for jj=ii+1:nCharges
            rj = rVec(:,jj);
            qj = q(jj);

            rij = ri - rj;
            dr = norm(rij);
            Usum = Usum +qi*qj/dr;
        end
    end
    Udirect(kk) = kB4pieps0_Inv* charge^2/angst* Usum;
end
plot(nList, Udirect, 'ob', 'linewidth', 2);

%% Ewald summation
% nReal = 10;
% nImag = 10;
% xi = 5.6/L*10; % Ewald parameter
% nCharges = 2;
% rVec = [-0.1 0.1;0 0;0 0];
% q=[-1;1];
% 
% USelf    = -kB4pieps0_Inv* charge^2/angst* sum(q.^2)*xi/sqrt(pi);
% UReal    =  kB4pieps0_Inv* charge^2/angst* realSum(a1,a2,a3,rVec,q,nCharges,nReal,L,xi);
% UFourier =  kB4pieps0_Inv* charge^2/angst* waveSum(a1,a2,a3,rVec,q,nCharges,nImag,L,xi);
% 
% Utot = UReal + UFourier + USelf;
%
% %% Miscellaneous
% Visualize charges in space
% for ii=1:nCharges
%     scatter3(rVec(1,ii), rVec(2,ii), rVec(3,ii));
%     hold on;
% end
toc