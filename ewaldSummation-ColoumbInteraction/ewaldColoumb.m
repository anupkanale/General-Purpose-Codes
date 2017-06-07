% Program to illustrate the use of Ewald summations to calculate total
% electrostatic energy of a system of charges due to Coloumbic interaction
clear; clc;
close all;

N = 6; % Number of particles per box
eps0 = 8.8e-12; % vacuum permitivitty
const = 4*pi*eps0*2;
L = 1; % unit cell
sigma = 0.4; % SD of superimposed gaussian cloud

q1 = ones(N/2,1); 
q2 = -ones(N/2,1); % satisfy neutral charge condition
q = cat(1,q1,q2);
nBoxes = 5^3; % number of periodic boxes
a1 = [1;0;0]; % lattice vectors
a2 = [0;1;0];
a3 = [0;0;1];

% Assign positions to charges
rng('default');
rng(1); % initialise seed to generate repeatable random number distribution
r = rand(3,N);

%% Direct Computation
UDirect = directSum(a1,a2,a3,r,q,N,nBoxes,L,eps0);

%% Ewald summation
URealSum    = eps0*realSum(a1,a2,a3,r,q,N,nBoxes,L,eps0,sigma);
UFourierSum = eps0*waveSum(a1,a2,a3,r,q,N,nBoxes,L,eps0,sigma);
USelf = eps0*selfInterac(q,eps0,sigma);
Utot =URealSum + UFourierSum - USelf;

err = (Utot-UDirect)/UDirect;
disp(USelf);