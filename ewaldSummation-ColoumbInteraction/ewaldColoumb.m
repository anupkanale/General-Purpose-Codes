% Program to illustrate the use of Ewald summations to calculate total
% electrostatic energy of a system of charges due to Coloumbic interaction
clear; %clc;

angst = 1e-10;
charge = 0.4238*1.6e-19;

%% read data from external file
filename = 'validationData/spce_sample_config_periodic1.txt';
data = readtable(filename);
data(:,1) = []; % delete serial number

r = transpose(table2array(data(:, 1:3)))*angst;
ions = char(table2array(data(:,4)));
q = zeros(length(ions),1);
for ii=1:length(ions)
    if ions(ii)=='O'
        q(ii) = -2;
    elseif ions(ii)=='H'
        q(ii) = 1;
    end
end
q = q*charge;

N = length(ions); % number of charged particles
M=N/3; % number of molecules
fHandle = fopen(filename, 'r');
firstLine = sscanf(fgetl(fHandle), '%f');
L = firstLine(1)*angst;
fclose(fHandle);

%% Problem setup-- Assign charges & positions, set parameters
nReal = 4;
nImag = 0;
alpha = 5.6/L; % Ewald parameter
a1 = [1;0;0]; % lattice vectors
a2 = [0;1;0];
a3 = [0;0;1];
eps0 = 8.85e-12;
kB = 1.38e-23;
const = 1/(4*pi*eps0*kB);

%% Ewald summation

USelf       = -const*  sum(q.^2)*alpha/sqrt(pi);
% UIntra      = -const* intraSum(r,q,alpha,M,L);
% URealSum    = const* realSum(a1,a2,a3,r,q,N,nReal,L,alpha);
% UFourierSum = const* waveSum(a1,a2,a3,r,q,N,nImag,L,alpha);

% Utot =URealSum + UFourierSum + USelf +UIntra;

nList = 0:5;
len = length(nList);
alphaList = 5.6/L*(1:0.5:5);
UPlot = zeros(length(nList),length(alphaList));

for jj=1:length(alphaList)
    alpha= alphaList(jj);
    parfor ii=1:len 
        nReal = nList(ii);
        URealSum  = const/2* realSum(a1,a2,a3,r,q,N,nReal,L,alpha);
        UPlot(ii,jj) = URealSum;
    end
end

%% Plot
close all;
set(figure(), 'position', [50 50 1000 800]);
legendInfo = cell(length(alphaList),1);
for jj=1:length(alphaList)
    legendInfo{length(alphaList)-jj+1} = strcat('\alpha = ', num2str(alphaList(jj)*1e-10));
    plot(nList, UPlot(:,jj), 'o-', 'linewidth', 1.5);
    hold on;
end
title('Convergence plot for real sum');
xlabel('Number of periodic shells'); ylabel('RealSum');
legend(legendInfo, 'location', 'best');

% %% Visualize primary cell
% %-----------------------
% set(figure(1), 'position', [500 500 1000 1000]);
% set(gcf,'color','w');
% for ii=1:300
%     if ions(ii)=='O'
%         scatter3(r(1,ii), r(2,ii), r(3,ii), 'filled', 'r');
%     elseif ions(ii)=='H'
%         scatter3(r(1,ii), r(2,ii), r(3,ii), 'filled', 'markerfacecolor', [0.5 0.5 0.5]);
%     end
%     hold on;
% %     text(r(1,ii),r(2,ii),r(3,ii),num2str(ii))
% end
% legend('Oxygen atom', 'Hydrogen atom', 'location', 'best')
% grid off;
% set(gca, 'Box', 'on', 'BoxStyle', 'full');
% set(gca,'xtick',[], 'ytick',[], 'ztick',[])