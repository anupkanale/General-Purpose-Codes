% Program to illustrate the use of Ewald summations to calculate total
% electrostatic energy of a system of charges due to Coloumbic interaction
clear; %clc;
tic

%% Data from file
angst = 1e-10;
charge = 0.4238*1.6e-19;

fileNum = 1;
filename = strcat('validationData/spce_sample_config_periodic', num2str(fileNum), '.txt');
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
 
% %% Ewald summation
% USelf    = -kB4pieps0_Inv* charge^2/angst* sum(q.^2)*alpha/sqrt(pi);
% UIntra   = -kB4pieps0_Inv* charge^2/angst* intraSum(r,q,alpha,M,L);
% UReal    =  kB4pieps0_Inv* charge^2/angst* realSum(a1,a2,a3,r,q,N,nReal,L,alpha);
% UFourier =  kB4pieps0_Inv* charge^2/angst* waveSum(a1,a2,a3,r,q,N,nImag,L,alpha);
% 
% Utot = UReal + UFourier + USelf + UIntra;
 
% %% Visualize primary cell
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
% 
% 
% %%
% set(figure(), 'position', [0 0 1000 800]);
% 
% subplot(2,2,1);
% UReal_ref = [-5.58889e+05	-1.19295e+06	-1.96297e+06	-3.57226e+06];
% plot(N,UReal/1e6, '-ro', N, UReal_ref/1e6, '-bo', 'linewidth', 1.5);
% legend('myValue', 'reference');
% title('Real space sum');
% xlabel('Number of charges'); ylabel('UReal');
% set(gca,'color', 'w');
% 
% 
% subplot(2,2,2);
% UFourier_ref = [6.27009e+03	6.03495e+03	5.24461e+03	7.58785e+03];
% plot(N,UFourier/1e3, '-ro', N,UFourier_ref/1e3, '-bo','linewidth', 1.5);
% legend('myValue', 'reference');
% title('Fourier space sum');
% xlabel('Number of charges'); ylabel('UFourier');
% set(gca, 'color', 'w');
% 
% 
% subplot(2,2,3);
% USelf_ref = [-2.84469e+06	-5.68938e+06	-8.53407e+06	-1.42235e+07];
% plot(N,USelf/1e6, '-ro', N,USelf_ref/1e6, '-bo','linewidth', 1.5);
% legend('myValue', 'reference');
% title('Sel-interac energy ');
% xlabel('Number of charges'); ylabel('USelf');
% set(gca,'color', 'w');
% 
% 
% subplot(2,2,4);
% UIntra_ref = [2.80999e+06	5.61998e+06	8.42998e+06	1.41483e+07];
% plot(N,UIntra/1e6, '-ro', N,UIntra_ref/1e6, '-bo','linewidth', 1.5);
% legend('myValue', 'reference');
% title('Intra-molecular energy');
% xlabel('Number of charges'); ylabel('UIntra');
% set(gca,'color', 'w');


%% Test for convergence of real sum
alphaList = 5.6/L*(1:6);

nList1 = 0:8;
UPlot1 = zeros(length(nList1),length(alphaList));

for jj=1:length(alphaList)
    alpha= alphaList(jj);
    parfor ii=1:length(nList1)
        nReal = nList1(ii);
        UReal  = kB4pieps0_Inv*charge^2/angst* realSum(a1,a2,a3,r,q,N,nReal,L,alpha);
        UPlot1(ii,jj) = UReal;
    end
end
%%
set(figure(), 'position', [50 50 1000 800]);
for jj=1:length(alphaList)
     plot(nList1, UPlot1(:,jj), 'o-', 'linewidth', 1.5);
     hold on;
end
title('Convergence plot for real space sum', 'fontsize', 20);
xlabel('Number of periodic shells', 'fontsize', 20); ylabel('realSum', 'fontsize', 20);
% legn = legend('\xi_1=0.28','\xi_2=0.56','\xi_3=0.84');
% legn.FontSize = 20;
% legn.Location = 'best';
set(gca,'fontsize',20, 'color', 'w')

%% Test for convergence of Imag sum
nList2 = 10:5:30;
len = length(nList2);
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
% legn = legend('\xi_1=0.28','\xi_2=0.56','\xi_3=0.84');
% legn.FontSize = 20;
% legn.Location = 'best';
set(gca,'fontsize',20, 'color', 'w')
 
%%
Utot = zeros(length(alphaList),1);
Uself = zeros(length(alphaList),1);
for jj=1:length(alphaList)
    Uself(jj) = kB4pieps0_Inv* charge^2/angst*sum(q.^2)*alphaList(jj)/sqrt(pi);
    Utot(jj) = UPlot1(5,jj) + UPlot2(5,jj) - Uself(jj);
end
 
set(figure(), 'position', [50 50 1000 800]);
plot(alphaList,-UPlot1(5,:), 'o-', alphaList,-UPlot2(5,:),'o-', ...
     alphaList, Uself,'o-',  alphaList, -Utot, 'o-', 'linewidth',2);
title('Total energy plot','fontsize',20);
xlabel('\xi', 'fontsize', 20); ylabel('U', 'fontsize', 20);
% legn = legend('U_{real}','U_{Fourier}','U_{self}', 'U_{tot}');
% legn.FontSize = 20;
% legn.Location = 'best';
set(gca,'fontsize',20, 'color', 'w')

save('dddummy')

toc
