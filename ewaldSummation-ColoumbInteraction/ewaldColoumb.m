% Program to illustrate the use of Ewald summations to calculate total
% electrostatic energy of a system of charges due to Coloumbic interaction
clear; clc;

%% read data from external file
filename = 'validationData/spce_sample_config_periodic1.txt';
data = readtable(filename);

data(:,1) = []; % delete serial number
r = transpose(table2array(data(:, 1:3)))*1e-10;

% scatter3(r(1,:)+10e-10, r(2,:)+10e-10, r(3,:)+10e-10);
% ylabel('Yaxis');
% zlabel('Z axis');

ions = char(table2array(data(:,4)));
q = zeros(length(ions),1);
for ii=1:length(ions)
    if ions(ii)=='O'
        q(ii) = -2*0.4238*1.6e-19;
    elseif ions(ii)=='H'
        q(ii) = 1*0.4238*1.6e-19;
    end
end

N = length(ions);
fHandle = fopen(filename, 'r');
firstLine = sscanf(fgetl(fHandle), '%f');
L = firstLine(1)* 1e-10;
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

USelf       = const* sum(q.^2)*alpha/sqrt(pi);
URealSum    = const* realSum(r,q,N,nReal,L,alpha);
% UFourierSum = const* waveSum(a1,a2,a3,r,q,N,nImag,L,alpha);

% Utot =URealSum + UFourierSum - USelf;

for ii=1:N
% nRealList = (3:2:35).^3;
% URealSum = zeros(length(nRealList),1);
% for ii=1:length(nRealList)
%     nReal = nRealList(ii);
%     URealSum(ii)    = realSum(a1,a2,a3,r,q,N,nReal,L,alpha);
% end
% set(figure(1), 'position', [3000 1000 800 700]);
% plot(nRealList, abs(URealSum), 'linewidth', 2);

% nImagList = (5).^3;
% UFourierSum = zeros(length(nImagList),1);
% for ii=1:length(nImagList)
%     nImag = nImagList(ii);
%     UFourierSum(ii) = waveSum(a1,a2,a3,r,q,N,nImag,L,alpha);
% end
% set(figure(2), 'position', [3000 50 800 700]);
% plot(nImagList, UFourierSum,'linewidth', 2);
end
