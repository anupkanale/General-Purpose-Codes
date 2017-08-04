%% Program to simulate flow past a stokeslet using direct computation
clear; %close all;

eta = 1;

% grid points
res = 15;
nX = res;
nY = res*2+1;
nZ = 1; %res;
nPoints = nX*nY*nZ;
rx = linspace(-1, 1, nX);
ry = linspace(-1, 1, nY);
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
nStokes = 2; % number of stokeslets
fVec = zeros(3,nStokes);
r_fVec = [-0.1 0.1 ; 0 0 ; 0 0];
fVec(:,1) = 5*[-1; 0; 0];
fVec(:,2) = 5*[1; 0; 0];
% fVec(:,3) = [-1; 0; 0];
% fVec(:,4) = [1; 0; 0];

%% compute flowfield
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
save('direct_Stokes_data');

clf(figure(1));
set(figure(1), 'position', [50 250 850 700])
for ii=1:nStokes
    plot(r_fVec(1,ii), r_fVec(2,ii), 'or','markersize', 10, 'markerfacecolor', 'r');
    hold on;
end

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
caxis([0 5]);

quiv = quiver(rx,ry,u3d',v3d');
quiv.Color = 'red';
quiv.LineWidth = 1.5;
% axis equal;
grid minor;
xlim([-1.01 1.01]); ylim([-1.01 1.01]); % zlim([-0.5*L(3), 0.5*L(3)])
% view(2);

dummy = velVec';
dummy2 = u3d';
dummy3 = v3d';
title('Direct Computation');

% 3D plotting
% [x,y,z] = meshgrid(rx,ry,rz);
% quiv = quiver3(x,y,z,u3d,v3d,w3d);
% quiv.Color = 'White';
% quiv.LineWidth = 1.0;
% Slice plot
% xslice = 1;
% yslice = 1;
% zslice = 0.5;
% slice(x,y,z,velMag3d,xslice,yslice,zslice)


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

