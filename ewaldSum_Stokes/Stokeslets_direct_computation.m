%% Program to simulate flow past a stokeslet using direct computation
clear; %close all;

eta = 1;

% grid points
res = 13;
nX = res;
nY = res;
nZ = res;
nPoints = nX*nY*nZ;
lim = 1;
rx = linspace(0,lim,nX);
ry = linspace(0,lim,nY);
rz = linspace(0,lim*1000,nZ);
rVec = zeros(3,nPoints);
for kk=1:nZ
    for ii=1:nX
        for jj=1:nY
            pointNum = (ii-1)*nX + jj + (kk-1)*nX*nY;
            rVec(:,pointNum) = [rx(jj); ry(ii); rz(kk)];
        end
    end
end

% assign point force location
nStokes = 2; % number of stokeslets
fVec = zeros(3,nPoints);
offset = 5;
fLoc = [(nPoints+1)/2 + offset, (nPoints+1)/2 - offset];
fVec(:,fLoc(1)) = [0.1; 0; 0];
fVec(:,fLoc(2)) = [-0.1; 0; 0];

%% compute flowfield
velVec = zeros(3,nPoints);
for ii=1:nPoints
    rVecI = rVec(:,ii);
    sum = 0;
    for jj=1:nStokes
        rVecJ = rVec(:,fLoc(jj));
        rij = rVecI - rVecJ;
        
        hiMat = calcOseenTensor(rij,eta);
        sum = sum + hiMat *fVec(:,fLoc(jj));
    end
    velVec(:,ii) = sum;
end

%% Post-processing
save('direct_Stokes_data');
plotGen("direct");


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
