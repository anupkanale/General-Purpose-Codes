% Program to use Ewald summations to calculate total electrostatic energy
% due to Coloumbic interaction
clear; clc;
close all;

N = 2; % Number of particles
eps0 = 8.8e-12; % vacuum permitivitty
const = 4*pi*eps0*2;
q1 = ones(N/2,1); % satisfy neutral charge condition
q2 = -ones(N/2,1);
q = cat(1,q1,q2);
sigma = 0.8; % SD of gaussian superimposed cloud

box = 3:2:25;
box = box';
Uplot = zeros(length(box),1);
for trial=1:length(box)
nBoxes = box(trial)^3; % number of periodic boxes
L = 1; % unit cell

% Assign positions to charges
rng('default');
rng(1); % initialise seed to generate repeatable random number distribution
r = rand(3,N);

%% Real Space
% generate periodic box

a1 = [1;0;0]; a2 = [0;1;0]; a3=[0;0;1]; % lattice vectors
[pVec] = makePeriodicBox(a1, a2, a3, L, nBoxes);

% Short-range/Real-space sum
URealSum = 0;
for kk=1:nBoxes % periodic box #
    for ii=1:N
    for jj=1:N
        if( norm(pVec(:,kk))==0 && jj==ii ), URealSum = URealSum + 0;
        else
            dist = norm( r(:,ii) - r(:,jj) + pVec(:,kk) );
            Uij = q(ii)*q(jj)/dist * erfc(dist/sqrt(2)/sigma);
            URealSum = URealSum + Uij;
        end
    end
    end
end
% URealSum= URealSum/const;

%% Fourier Space

%Reciprocal Lattice vectors
Vol = dot(a1,cross(a2,a3));
b1 = cross(a2,a3)/Vol;
b2 = cross(a3,a1)/Vol;
b3 = cross(a1,a2)/Vol;

[pKVec] = makePeriodicBox(b1, b2, b3, L, nBoxes);
% Fourier Space summation
UFourierSum = 0;
for kk=1:nBoxes    
    kVec = pKVec(:,kk);
    kMag = norm(kVec);
    if kMag~=0
        for ii=1:N
        for jj=1:N
            rij = r(:,ii) - r(:,jj);
            UijF = q(ii)*q(jj) /kMag^2 * exp(-1j* dot(kVec,rij) );
            UijF = UijF * exp(-sigma^2 * kMag^2/2);
            UFourierSum = UFourierSum + UijF;
        end
        end
    end
end
UFourierSum = UFourierSum/(2*Vol*eps0);

%% Self Interaction
Uself=0;
for ii=1:N
    Uself = Uself + q(ii)^2;
end
Uself = Uself/(4*pi*eps0*sigma*sqrt(2*pi));

%% Total Energy
Utot = URealSum + UFourierSum - Uself;
Uplot(trial) = URealSum;
end

close all;
err = zeros(length(box),1);
for trial=2:(length(box))
    err(trial) = abs(Uplot(trial)-Uplot(trial-1));
end
plot(box(2:end), err(2:end), 'o-', 'LineWidth', 1.2);
xlabel('N_b','FontSize', 10, 'FontWeight', 'bold');
ylabel('Error','FontSize', 10, 'FontWeight', 'bold');