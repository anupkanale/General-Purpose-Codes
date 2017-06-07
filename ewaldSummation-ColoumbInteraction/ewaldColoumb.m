% Program to illustrate the use of Ewald summations to calculate total
% electrostatic energy of a system of charges due to Coloumbic interaction
clear; clc;
close all;

N = 6; % Number of particles per box
alpha = 0.1; % Ewald parameter
L = 1; % unit cell

q1 = ones(N/2,1);
q2 = -ones(N/2,1); % satisfy neutral charge condition
q = cat(1,q1,q2);
nBoxes = 25^3; % number of periodic boxes
a1 = [1;0;0]; % lattice vectors
a2 = [0;1;0];
a3 = [0;0;1];

% Assign positions to charges
rng('default');
rng(1); % initialise seed to generate repeatable random number distribution
r = rand(3,N);

%% Direct Computation
UDirect = directSum(a1,a2,a3,r,q,N,nBoxes,L);

%% Ewald summation
URealSum    = realSum(a1,a2,a3,r,q,N,nBoxes,L,alpha);
UFourierSum = waveSum(a1,a2,a3,r,q,N,nBoxes,L,alpha);
USelf = selfInterac(q,alpha);
Utot =URealSum + UFourierSum - USelf;

err = (Utot-UDirect)/UDirect;