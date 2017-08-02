%% Program to simulate flow past a stokeslet using direct computation
clear; %close all;

eta = 1;

% grid points
res = 9;
nX = res;
nY = res;
nZ = 1; %res;
nPoints = nX*nY*nZ;
L = [20000; 20000; 20000];
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
nStokes = 4; % number of stokeslets
fVec = zeros(3,nStokes);
r_fVec = [-0.1 0.1 0 0; 0 0 -0.1 0.1; 0 0 0 0];
fVec(:,1) = [0; -1; 0];
fVec(:,2) = [0; 1; 0];
fVec(:,3) = [-1; 0; 0];
fVec(:,4) = [1; 0; 0];

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
set(figure(1), 'position', [1350 550 550 400])
for ii=1:nStokes
    scatter3(r_fVec(1,ii), r_fVec(2,ii), r_fVec(3,ii), 'filled'); % 'or','markersize', 10, 'markerfacecolor', 'r');
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

% pcolor(rx,ry,velMag2d);
% shading interp;
% colorbar;
scale = 500;
quiv = quiver(rx,ry,u2d*scale,v2d*scale);
quiv.Color = 'blue';
quiv.LineWidth = 1.5;
% axis equal;
grid minor;
xlim([-1 1]); ylim([-1 1]); % zlim([-0.5*L(3), 0.5*L(3)])
view(2);

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

