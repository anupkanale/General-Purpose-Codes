% Program to illustrate the use of Ewald summations to calculate total
% electrostatic energy of a system of charges due to Coloumbic interaction
clear; %clc;

angst = 1e-10;
charge = 0.4238*1.6e-19;

%% read data from external file
filename = 'validationData/spce_sample_config_periodic1.txt';
data = readtable(filename);
data(:,1) = []; % delete serial number

r = transpose(table2array(data(:,1:3)));
ions = char(table2array(data(:,4)));
q = zeros(length(ions),1);
for ii=1:length(ions)
    if ions(ii)=='O'
        q(ii) = -2;
    elseif ions(ii)=='H'
        q(ii) = 1;
    end
end

N = length(ions); % number of charged particles
M=N/3; % number of molecules
fHandle = fopen(filename, 'r');
firstLine = sscanf(fgetl(fHandle), '%f');
L = firstLine(1);
fclose(fHandle);

%% Problem setup-- Assign charges & positions, set parameters
nReal = 0;
nImag = 5;
alpha = 5.6/L; % Ewald parameter
a1 = [1;0;0];  % lattice vectors
a2 = [0;1;0];
a3 = [0;0;1];
eps0 = 8.85e-12;
kB = 1.38e-23;
kB4pieps0_Inv = 1/(kB*4*pi*eps0);

%% Ewald summation
% USelf    = -kB4pieps0_Inv* charge^2/angst* sum(q.^2)*alpha/sqrt(pi);
% UIntra   = -kB4pieps0_Inv* charge^2/angst* intraSum(r,q,alpha,M,L);
% UReal    =  kB4pieps0_Inv* charge^2/angst* realSum(a1,a2,a3,r,q,N,nReal,L,alpha);
% UFourier =  kB4pieps0_Inv* charge^2/angst* waveSum(a1,a2,a3,r,q,N,nImag,L,alpha);
% 
% Utot = UReal + UFourier + USelf + UIntra;

%% Visualize primary cell
% %-----------------------
% set(figure(1), 'position', [500 500 1000 1000]);
% set(gcf,'color','w');
% for ii=1:300
%     if ions(ii)=='O'
%         scatter3(r(1,ii), r(2,ii), r(3,ii), 'filled', 'r');
%     elseif ions(ii)=='H'
%         scatter3(r(1,ii), r(2,ii), r(3,ii), 'filled',...
%             'markerfacecolor', [0.5 0.5 0.5]);
%     end
%     hold on;
% %     text(r(1,ii),r(2,ii),r(3,ii),num2str(ii))
% end
% legend('Oxygen atom', 'Hydrogen atom', 'location', 'best')
% grid off;
% set(gca, 'Box', 'on', 'BoxStyle', 'full');
% set(gca,'xtick',[], 'ytick',[], 'ztick',[])

%% Test for convergence of real sum
nList = [0,1,2,4,8];
len = length(nList);
alphaList = 5.6/L*(1:3);
UPlot1 = zeros(length(nList),length(alphaList));

for jj=1:length(alphaList)
    alpha= alphaList(jj);
    parfor ii=1:len
        nReal = nList(ii);
        UReal  = kB4pieps0_Inv*charge^2/angst* realSum(a1,a2,a3,r,q,N,nReal,L,alpha);
        UPlot1(ii,jj) = UReal;
    end
end

%%
close all;
set(figure(), 'position', [50 50 1000 800]);
for jj=1:length(alphaList)
    plot(nList, UPlot1(:,jj), 'o-', 'linewidth', 1.5);
    hold on;
end
title('Convergence plot for real space sum', 'fontsize', 20);
xlabel('Number of periodic shells', 'fontsize', 20); ylabel('realSum', 'fontsize', 20);
legn = legend('\xi_1=0.28','\xi_2=0.56','\xi_3=0.84');
legn.FontSize = 20;
legn.Location = 'best';
set(gca,'fontsize',20)

%% Test for convergence of Imag sum
nList2 = 4:2:16;
len = length(nList2);
alphaList = 5.6/L*(1:3);
UPlot2 = zeros(length(nList2),length(alphaList));

for jj=1:length(alphaList)
    alpha= alphaList(jj);
    parfor ii=1:len
        nImag = nList2(ii);
        UFourier  = kB4pieps0_Inv*charge^2/angst* waveSum(a1,a2,a3,r,q,N,nImag,L,alpha);
        UPlot2(ii,jj) = UFourier;
    end
end

%%
set(figure(), 'position', [50 50 1000 800]);
for jj=1:length(alphaList)
    plot(nList2, UPlot2(:,jj), 'o-', 'linewidth', 1.5);
    hold on;
end
title('Convergence plot for k-space sum', 'fontsize', 20);
xlabel('Number of periodic shells', 'fontsize', 20); ylabel('waveSum', 'fontsize', 20);
legn = legend('\xi_1=0.28','\xi_2=0.56','\xi_3=0.84');
legn.FontSize = 20;
legn.Location = 'best';
set(gca,'fontsize',20)

%%
save('ewald_coloumb_data')

%%
for ii=1:3
    Uself(ii) = kB4pieps0_Inv* charge^2/angst*sum(q.^2)*alphaList(ii)/sqrt(pi);
    Utot(ii) = UPlot1(5,ii) + UPlot2(7,ii) - Uself(ii);
end
% plot(alphaList,Uself, 'linewidth',1.5);
% ylim([-2e7 -1e7])


plot(alphaList,-UPlot1(5,:), alphaList,-UPlot2(7,:), alphaList,Uself, alphaList, -Utot, 'linewidth',2);
legn = legend('U_{real}','U_{Fourier}','U_{self}', 'U_{tot}');
legn.FontSize = 20;
legn.Location = 'best';
set(gca,'fontsize',20)
xlabel('\xi', 'fontsize', 20); ylabel('U', 'fontsize', 20);