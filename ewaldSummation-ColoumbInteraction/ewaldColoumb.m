% Program to illustrate the use of Ewald summations to calculate total
% electrostatic energy of a system of charges due to Coloumbic interaction
clear; clc;

N = 4; % Number of particles per box
alpha = 0.8; % Ewald parameter
L = 1; % unit cell

q1 = ones(N/2,1);
q2 = -ones(N/2,1); % satisfy neutral charge condition
q = cat(1,q1,q2);
% nBoxes = 25^3; % number of periodic boxes
nReal = 25^3;
nImag = 25^3;
a1 = [1;0;0]; % lattice vectors
a2 = [0;1;0];
a3 = [0;0;1];

% Assign positions to charges
rng('default');
rng(1); % initialise seed to generate repeatable random number distribution
r = rand(3,N);

%% Direct Computation
% UDirect = directSum(a1,a2,a3,r,q,N,nBoxes,L);

%% Ewald summation

% nRealList = (3:2:35).^3;
% URealSum = zeros(length(nRealList),1);
% for ii=1:length(nRealList)
%     nReal = nRealList(ii);
%     URealSum(ii)    = realSum(a1,a2,a3,r,q,N,nReal,L,alpha);
% end
% plot(nRealList, abs(URealSum), 'linewidth', 2);

nImagList = (3:2:15).^3;
UFourierSum = zeros(length(nImagList),1);
for ii=1:length(nImagList)
    nImag = nImagList(ii);
    UFourierSum(ii) = waveSum(a1,a2,a3,r,q,N,nImag,L,alpha);
end
plot(nImagList, UFourierSum,'linewidth', 2);

% USelf = selfInterac(q,alpha);
% Utot =URealSum + UFourierSum - USelf;