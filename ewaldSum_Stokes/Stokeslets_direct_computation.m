clear;
close all;

eta = 1;
% N = 2; % number of stokeslets

% grid points
nCols = 11;
nRows = 11;
nPoints = nCols*nRows;

lim = 1;
rx = linspace(-lim,lim,nCols);
ry = linspace(-lim,lim,nRows);
rVec = zeros(2,nPoints);
for ii=1:nRows
    for jj=1:nCols
        pointNum = (ii-1)*nCols+jj;
        rVec(:,pointNum) = [rx(jj); ry(ii)];
    end
end

% point force location
fVec = zeros(2,nPoints);
% fVec(:,(nPoints+1)/2) = [1; 0];
fLoc = [(nPoints+1)/2 + 2, (nPoints+1)/2 - 2];

fVec(:,fLoc(1)) = [1; 0];
fVec(:,fLoc(2)) = [-1; 0];

% Compute flowfield
uVec = zeros(2,nPoints);
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
    uVec(:,ii) = sum;
end

%% Plotting
%----------------
scale = 500;
set(figure, 'Position', [2700, 1000, 1000, 700]);
title('Streamlines for slow past a Stokeslet');
for ii=1:length(fLoc)
    plot(rVec(1,fLoc(ii)), rVec(2,fLoc(ii)), 'markersize', 5);
    hold on;
end

% Convert to format appropriate for plotting
for pointNum=1:nPoints
    ii = floor((pointNum-1)/nCols) + 1;
    jj = rem(pointNum-1,nCols) + 1;
    u(ii,jj) = uVec(1,pointNum);
    v(ii,jj) = uVec(2,pointNum);
end

% Plot flowfield
q = quiver(rx,ry,u*scale,v*scale);
q.Color = 'blue';
q.LineWidth = 1.5;

% starty = linspace(-0.9,1.1,nCols);
% startx = -1*ones(size(starty));
% h = streamline(rx,ry,u,v,startx,starty);
% set(h, 'LineWidth', 2, 'Color', 'red');
axis equal;
grid minor;

%% Verification Lines
dummy = rVec';
dummy2 = fVec';
% % Check if the point numbers are being created correctly
% for pointNum=1:nPoints
%     ii = floor((pointNum-1)/nCols) + 1;
%     jj = rem(pointNum-1,nCols) + 1;
%     pointLocation(ii,jj) = (ii-1)*nCols+jj;
% end