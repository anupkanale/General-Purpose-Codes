% Program to illustrate the use of Ewald summations to calculate total
% electrostatic energy of a system of charges due to Coloumbic interaction
clear; %clc;
close all;

angst = 1e-10;
charge = 0.4238*1.6e-19;

%% read data from external file
filename = 'validationData/spce_sample_config_periodic1.txt';
data = readtable(filename);

data(:,1) = []; % delete serial number
r = transpose(table2array(data(:, 1:3)))*angst;

% set(figure(1), 'position', [2500 0 1000 1000]);
% % scatter3(r(1,:), r(2,:), r(3,:));
% for i=69:75%1:300
% scatter3(r(1,i), r(2,i), r(3,i));
% hold on;
% text(r(1,i),r(2,i),r(3,i),num2str(i))
% end
% ylabel('Yaxis');
% zlabel('Z axis');

ions = char(table2array(data(:,4)));
q = zeros(length(ions),1);
for ii=1:length(ions)
    if ions(ii)=='O'
        q(ii) = -2*charge;
    elseif ions(ii)=='H'
        q(ii) = 1*charge;
    end
end

N = length(ions); % number of charged particles
M=N/3; % number of molecules
fHandle = fopen(filename, 'r');
firstLine = sscanf(fgetl(fHandle), '%f');
L = firstLine(1)*angst;
fclose(fHandle);

%% Problem setup-- Assign charges & positions, set parameters
nReal = 1; % number of terms in real sum
nImag = 5^3; % number of terms in fourier sum
alpha = 5.6/L; % Ewald parameter
a1 = [1;0;0]; % lattice vectors
a2 = [0;1;0];
a3 = [0;0;1];
eps0 = 8.85e-12;
kB = 1.38e-23;
const = 1/(4*pi*eps0*kB);

%% Ewald summation

USelf       = -const*  sum(q.^2)*alpha/sqrt(pi);
UIntra = -const* intraSum(r,q,alpha,M,L);
% URealSum    = const* realSum(r,q,N,nReal,L,alpha);
UFourierSum = const* waveSum(a1,a2,a3,r,q,N,nImag,L,alpha);

% Utot =URealSum + UFourierSum + USelf +UIntra;

