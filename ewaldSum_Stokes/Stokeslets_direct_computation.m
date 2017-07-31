clear; close all;

eta = 1;
% N = 2; % number of stokeslets

% grid points
res = 7;
nX = res;
nY = res;
nZ = res;
nPoints = nX*nY*nZ;
lim = 1;
rx = linspace(-lim,lim,nX);
ry = linspace(-lim,lim,nY);
rz = linspace(-lim,lim,nZ);
rVec = zeros(3,nPoints);
for kk=1:nZ
    for ii=1:nX
        for jj=1:nY
            pointNum = (ii-1)*nX + jj + (kk-1)*nX*nY;
            rVec(:,pointNum) = [rx(jj); ry(ii); rz(kk)];
        end
    end
end


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
% 
% 
% 
% point force location
fVec = zeros(3,nPoints);
fLoc = [(nPoints+1)/2 + 1, (nPoints+1)/2 - 1];

fVec(:,fLoc(1)) = [1; 0; 0];
fVec(:,fLoc(2)) = [-1; 0; 0];

% Compute flowfield
velVec = zeros(3,nPoints);
for ii=1:nPoints
    rVecI = rVec(:,ii);
    sum = 0;
    for jj=1:nPoints
        rVecJ = rVec(:,jj);
        rij = rVecI - rVecJ;
        
        if norm(rij)~=0
            hiMat = calcOseenTensor(rij,eta);
            sum = sum + hiMat *fVec(:,jj);
        end
    end
    velVec(:,ii) = sum;
end

%% Plotting
%----------------
scale = 500;
set(figure, 'Position', [2700, 1000, 1000, 700]);
title('Streamlines for slow past a Stokeslet');

[X,Y,Z] = sphere(4);
for ii=1:length(fLoc)
    scatter3(rVec(1,fLoc(ii)), rVec(2,fLoc(ii)), rVec(3,fLoc(ii)), 'filled');
    hold on;
end
xlim([-1 1]); ylim([-1 1]); zlim([-1 1])

% Convert to format appropriate for plotting
u3d = zeros(nX,nY,nZ);
v3d = zeros(nX,nY,nZ);
w3d = zeros(nX,nY,nZ);
for pointNum=1:nPoints
    kk = floor((pointNum-1)/(nX*nY)) + 1;
    jj = rem(pointNum-1,nX) + 1;
    ii = floor((pointNum-(kk-1)*nX*nY-1)/nX) + 1;
    u3d(ii,jj,kk) = velVec(1,pointNum);
    v3d(ii,jj,kk) = velVec(2,pointNum);
    w3d(ii,jj,kk) = velVec(3,pointNum);
end

% Plot flowfield in 2D
u2d = squeeze(u3d(:,:,(nZ+1)/2));
v2d = squeeze(v3d(:,:,(nZ+1)/2));
quiv = quiver(rx,ry,u2d*scale,v2d*scale);
quiv.Color = 'blue';
quiv.LineWidth = 1.5;
view(2);
axis equal;
grid minor;