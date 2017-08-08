%% Program to simulate flow past a stokeslet using direct computation
clear; %close all;

eta = 1/(8*pi);

% grid points
res = 101;
nX = res;
nY = res;
nZ = 1; %res;
nPoints = nX*nY*nZ;
L = 1;
box = L*[1; 1; 1];
rsize = L/2;
rx = linspace(-rsize,rsize,nX);
ry = linspace(-rsize,rsize,nY);
rz = 0;
rVec = zeros(3,nPoints);
for kk=1:nZ
    for jj=1:nY
        for ii=1:nX
            pointNum = ii + (jj-1)*nX + (kk-1)*nX*nY;
            rVec(:,pointNum) = [rx(ii); ry(jj); rz(kk)];
        end
    end
end

% point force location
% Layers=0 corresponds to non-periodic case.
layers = 1;
nStokes = 2*(2*layers+1)^3;
a1 = [1;0;0];
a2 = [0;1;0];
a3 = [0;0;1];
pMat = getP(a1,a2,a3,layers,box);

r_fVec = zeros(3,nStokes);
fVec = zeros(3,nStokes);

for ii=1:(2*layers+1)^3
    pVec = pMat(:,ii);
    r_fVec(:,2*ii-1) = [-0.1;0;0] + pVec;
    fVec(:,2*ii-1) = [-1; 0; 0];
    
    r_fVec(:,2*ii) = [0.1;0;0] + pVec;
    fVec(:,2*ii) = [1; 0; 0];
end

%% Compute flowfield
velVec = zeros(3,nPoints);
for ii=1:nPoints
    rVecI = rVec(:,ii);
    sum = 0;
    for jj=1:nStokes
        rVecJ = r_fVec(:,jj);
        rij = rVecI - rVecJ;
        
        hiMat = calcOseenTensor(rij,eta);
        sum = sum + hiMat *fVec(:,jj);
    end
    velVec(:,ii) = sum;
end

%% Post-processing
% save('direct_Stokes_data');

set(figure(), 'position', [150 0 1000 800], 'color', 'w')

% Convert to format appropriate for plotting
u3d = zeros(nX,nY,nZ);
v3d = zeros(nX,nY,nZ);
w3d = zeros(nX,nY,nZ);
velMag3d = zeros(nX,nY,nZ);
for pointNum=1:nPoints
    kk = floor((pointNum-1)/(nX*nY)) + 1;
    jj = floor((pointNum-(kk-1)*nX*nY-1)/nX) + 1;
    ii = rem(pointNum-1,nX) + 1;
    
    u3d(ii,jj,kk) = velVec(1,pointNum);
    v3d(ii,jj,kk) = velVec(2,pointNum);
    w3d(ii,jj,kk) = velVec(3,pointNum);
    velMag3d(ii,jj,kk) = norm(velVec(:,pointNum));
end

% Plot flowfield in 2D
u2d = squeeze(u3d(:,:,(nZ+1)/2));
v2d = squeeze(v3d(:,:,(nZ+1)/2));
velMag2d = squeeze(velMag3d(:,:,(nZ+1)/2));

[x,y] = meshgrid(rx,ry);
pcolor(x,y,velMag2d');
shading interp;
colorbar;
caxis([0 10]);

hold on;
quiv = quiver(rx,ry,u3d',v3d');
quiv.Color = 'red';
quiv.LineWidth = 1.5;
axis equal;
grid minor;
xlim([-rsize rsize]); ylim([-rsize rsize]); % zlim([-0.5*L(3), 0.5*L(3)])
% view(2);

hold on;
dummy = real(velVec');
dummy2 = rVec';
dummy3 = u3d';
dummy4 = v3d';
title('Direct Computation');

%%
for ii=1:nStokes
    plot(r_fVec(1,ii), r_fVec(2,ii), 'or','markersize', 10, 'markerfacecolor', 'r');
%     scatter3(r_fVec(1,ii), r_fVec(2,ii), r_fVec(3,ii));
%     hold on;
end

%% Verification
%-----------------
% % 3D plot to check if point numbers have been assigned correctly
% for kk=1:nZ
%     for ii=1:nX
%         for jj=1:nY
%             pointNum = (ii-1)*nX + jj + (kk-1)*nX*nY;
%             rVec(:,pointNum) = [rx(jj); ry(ii); rz(kk)];
%             
%             scatter3(rVec(1,pointNum), rVec(2,pointNum), rVec(3,pointNum),...
%                 'filled', 'markerfacecolor', [0.3*(kk-1) 0.3*(kk-1) 0.3*(kk-1)]);
%             string = strcat(num2str(ii), num2str(jj), num2str(kk), ':', num2str(pointNum));
%             text(rVec(1,pointNum), rVec(2,pointNum), rVec(3,pointNum), string);
%             hold on;
%         end
%     end
% end
% xlabel('X-axis'); ylabel('Y-axis'); zlabel('Z-axis');
% 
% % get i,j,k back from pointNum
% for pointNum=1:nPoints
%     kk = floor((pointNum-1)/(nX*nY)) + 1;
%     jj = rem(pointNum-1,nX) + 1;
%     ii = floor((pointNum-(kk-1)*nX*nY-1)/nX) + 1;
% end