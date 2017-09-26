%%
nReal = 3;
nImag = 9;
xiList = (2:20).^2;
a1 = [0; 0.5; 0.5];
a2 = [0.5; 0; 0.5];
a3 = [0.5; 0.5; 0];

% aa = [0,0.5];
% r = zeros(3,8);
% index = 1;
% for ii=1:2
%     for jj=1:2
%         for kk=1:2
%             r(:,index) = [aa(ii); aa(jj); aa(kk)];
%             index = index+1;
%         end
%     end
% end
% q = [1;-1;-1;1;-1;1;1;-1];
r = [0 0.25;0 0.25; 0 0.25];
q = [1 -1];
N = length(q);
L = 1;
loc = 1;

% parfor ii=1:length(xiList)
    xi = 7.07;
    phiReal = realSum(a1,a2,a3,r,q,N,nReal,L,xi,loc);
    phiFourier = waveSum(a1,a2,a3,r,q,N,nImag,L,xi,loc);
    phiSelf = 2*selfSum(q,N,xi,loc);

    phiTot = phiReal + phiFourier + phiSelf;
    madelung = real(phiTot/4*sqrt(3));
% end
%%
% figure()
% plot(xiList,-madelung, 'o-', 'linewidth', 2);