%% Program to simulate flow past a stokeslet using Ewald Summations
clear; %close all
tic

% grid points
% res = 5;
% nX = res;
% nY = res;
% nZ = 1;
% nPoints = nX*nY*nZ;
% L = 1;
% box = L*[1; 1; 1];
% rsize = L/2;
% rx = linspace(-rsize,rsize,nX);
% ry = linspace(-rsize,rsize,nY);
% rz = 0;
% rVec = zeros(3,nPoints);
% for kk=1:nZ
%     for jj=1:nY
%         for ii=1:nX
%             pointNum = ii + (jj-1)*nX + (kk-1)*nX*nY;
% %             rVec(:,pointNum) = [rx(ii); ry(jj); rz(kk)];
%         end
%     end
% end

nPoints = 1;
rVec = [0; 0; 0];
L = 1;
box = L*[1; 1; 1];
rsize = L/2;


% point force location
nStokes = 1; % number of stokeslets
fVec = zeros(3,nPoints);
r_fVec = [0.5; 0.5; 0.5];
fVec = [0.1;0.2;0.3];
% fVec(:,2) = [0.1;0.2;0.3];

% ewald summation parameters
a1 = [1;0;0];
a2 = [0;1;0];
a3 = [0;0;1];
nReal = 2;
nImag = 20;
xi = 15;

%% compute flowfield
velVec = zeros(3,nPoints);
uReal  = zeros(3,nPoints);
uFourier = zeros(3,nPoints);
uSelf = zeros(3,nPoints);
for mm=1:nPoints
    rVecM = rVec(:,mm);
    
%     uReal(:,mm) = realSum(nStokes,rVecM,r_fVec,fVec,xi,box,a1,a2,a3,nReal);
    uFourier(:,mm) = fourierSum(nStokes,rVecM,r_fVec,fVec,xi,box,a1,a2,a3,nImag);
%     uSelf(:,mm) = 4*xi/sqrt(pi) * fVec(:,);
    
    velVec(:,mm) = uReal(:,mm) + uFourier(:,mm) - uSelf(:,mm);
end


% %% Post-processing
% % save('ewald_Stokes_data');
% 
% % Convert to format appropriate for plotting
% u3d = zeros(nX,nY,nZ);
% v3d = zeros(nX,nY,nZ);
% w3d = zeros(nX,nY,nZ);
% velMag3d = zeros(nX,nY,nZ);
% for pointNum=1:nPoints
%     kk = floor((pointNum-1)/(nX*nY)) + 1;
%     jj = floor((pointNum-(kk-1)*nX*nY-1)/nX) + 1;
%     ii = rem(pointNum-1,nX) + 1;
%     
%     u3d(ii,jj,kk) = velVec(1,pointNum);
%     v3d(ii,jj,kk) = velVec(2,pointNum);
%     w3d(ii,jj,kk) = velVec(3,pointNum);
%     velMag3d(ii,jj,kk) = norm(velVec(:,pointNum));
% end
% 
% % Plot flowfield in 2D
% u2d = squeeze(u3d(:,:,(nZ+1)/2));
% v2d = squeeze(v3d(:,:,(nZ+1)/2));
% velMag2d = squeeze(velMag3d(:,:,(nZ+1)/2));
% 
% %% Plots
% set(figure(), 'position', [1050 1050 850 700]);
% 
% [x,y] = meshgrid(rx,ry);
% pcolor(x,y,velMag2d');
% shading interp;
% colorbar;
% caxis([0 5000]);
% 
% hold on;
% quiv = quiver(rx,ry,u3d',v3d');
% quiv.Color = 'red';
% quiv.LineWidth = 1.5;
% axis equal;
% grid minor;
% xlim([-rsize rsize]); ylim([-rsize rsize]); % zlim([-0.5*L(3), 0.5*L(3)])
% % view(2);
% title('Ewald method');
% 
% hold on;
% for ii=1:nStokes
%     plot(r_fVec(1,ii), r_fVec(2,ii), 'or','markersize', 10, 'markerfacecolor', 'r');
% end
% 
% dummyvel = real(velVec');
% dummyr = rVec';
% dummyu = u3d';
% dummyv = v3d';
% dummyreal = uReal';
% dummyfourier = uFourier';

toc