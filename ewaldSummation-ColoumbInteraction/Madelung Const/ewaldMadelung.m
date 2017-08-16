%%
nReal = 3;
nImag = 9;
xi = 25;
a1 = [1;0;0];
a2 = [0;1;0];
a3 = [0;0;1];

aa = [0,0.5];
r = zeros(3,8);
index = 1;
for ii=1:2
    for jj=1:2
        for kk=1:2
            r(:,index) = [aa(ii); aa(jj); aa(kk)];
            index = index+1;
        end
    end
end
q = [1;-1;-1;1;-1;1;1;-1];
N = length(q);
L = 1;
loc = 1;

phiReal = realSum(a1,a2,a3,r,q,N,nReal,L,xi,loc);
phiFourier = waveSum(a1,a2,a3,r,q,N,nImag,L,xi,loc);
phiSelf = 2*selfSum(q,N,xi,loc);

phiTot = phiReal + phiFourier + phiSelf;
madelung = phiTot/2;